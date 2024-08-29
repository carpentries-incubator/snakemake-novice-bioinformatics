###
# Snakefile you should have after completing episodes 01 to 07,
# but using wrappers from https://snakemake-wrappers.readthedocs.io
#
# To run a full MultiQC report on all samples, use:
#
# $ snakemake -j1 -p multiqc
###

# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# "rule all_counts" is no longer being used and has been removed to reduce clutter

# Generic read counter rule using wildcards and placeholders,
# which can count trimmed and untrimmed reads.
rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# Trim any FASTQ reads for base quality
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 22 -l 100 -o {output} <{input}"

# Kallisto quantification of one sample
# Modified to declare the whole directory as the output, and to capture all output to
# a log file.
rule kallisto_quant:
    output: directory("kallisto.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fastq = [ "trimmed/{sample}_1.fq", "trimmed/{sample}_2.fq" ],
    log:
        "logs/kallisto_quant_{sample}.log",
    wrapper:
        "v2.2.1/bio/kallisto/quant"

rule kallisto_index:
    output:
        index = "{strain}.kallisto_index",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    log:
        "{strain}.kallisto_log",
    wrapper:
        "v2.2.1/bio/kallisto/index"

rule fastqc:
    output:
        html = "{indir}.{myfile}_fastqc.html",
        zip  = "{indir}.{myfile}_fastqc.zip"
    input:  "{indir}/{myfile}.fq"
    log:
        "logs/fastqc/{indir}.{myfile}.log"
    wrapper:
        "v4.2.0/bio/fastqc"

rule salmon_quant:
    output:
        dir   = directory("salmon.{sample}"),
        quant = "salmon.{sample}/quant.sf",
        lib   = "salmon.{sample}/lib_format_counts.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
        r1   = "trimmed/{sample}_1.fq",
        r2   = "trimmed/{sample}_2.fq",
    log:
        "logs/salmon/{sample}.log",
    params:
        libtype = "A"
    wrapper:
        "v4.2.0/bio/salmon/quant"

rule salmon_index:
    output:
        idx = directory("{strain}.salmon_index")
    input:
        sequences = "transcriptome/{strain}.cdna.all.fa.gz"
    log:
        "logs/salmon/{strain}_index.log",
    wrapper:
        "v4.2.0/bio/salmon/index"

# A version of the MultiQC rule that ensures nothing unexpected is hoovered up by multiqc,
# by linking the files into a temporary directory.
# Note that this requires the *kallisto_quant* rule to be amended as above so that it has
# a directory as the output, with that directory containing the console log.
rule multiqc:
    output:
        "multiqc_out/multiqc.html",
        directory("multiqc_out/multiqc_data"),
    input:
        salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
    log:
        "logs/multiqc.log"
    wrapper:
        "v4.2.0/bio/multiqc"
