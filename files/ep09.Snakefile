###
# Snakefile you should have after completing episode 09, assuming you start with ep07.Snakefile
###

# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

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

# Kallisto quantification of one sample.
# Modified to declare the whole directory as the output.
rule kallisto_quant:
    output: directory("kallisto.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    threads: 4
    shell:
        """mkdir {output}
           kallisto quant -t {threads} -i {input.index} -o {output} {input.fq1} {input.fq2} >& {output}/kallisto_quant.log
        """

rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
        log = "{strain}.kallisto_log",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {output.log}"

rule fastqc:
    output:
        html = "{indir}.{myfile}_fastqc.html",
        zip  = "{indir}.{myfile}_fastqc.zip"
    input:  "{indir}/{myfile}.fq"
    shell:
        """fastqc -o . {input}
           mv {wildcards.myfile}_fastqc.html {output.html}
           mv {wildcards.myfile}_fastqc.zip  {output.zip}
        """

rule salmon_quant:
    output: directory("salmon.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    threads: 4
    shell:
        "salmon quant -p {threads} -i {input.index} -l A -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output}"

rule salmon_index:
    output:
        idx = directory("{strain}.salmon_index")
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "salmon index -t {input.fasta} -i {output.idx} -k 31"

# A version of the MultiQC rule that ensures nothing unexpected is hoovered up by multiqc,
# by linking the files into a temporary directory.
# Note that this requires the *kallisto_quant* rule to be amended as above so that it has
# a directory as the output, with that directory containing the console log.
rule multiqc:
    output:
        mqc_out = directory('multiqc_out'),
        mqc_in  = directory('multiqc_in'),
    input:
        salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
    shell:
        """mkdir {output.mqc_in}
           ln -snr -t {output.mqc_in} {input}
           multiqc {output.mqc_in} -o {output.mqc_out}
        """
