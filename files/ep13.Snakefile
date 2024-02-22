###
# Snakefile you should have after completing episode 13, assuming you start with ep07.Snakefile
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
        "echo $(( $(wc -l <{input:q}) / 4 )) > {output:q}"

LEN_READS_CMD = "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}"
rule lenreads:
    output: "{indir}.{myfile}.fq.len"
    input:  "{indir}/{myfile}.fq"
    shell:
        "awk {LEN_READS_CMD:q} {input:q} > {output:q}"


# Trim any FASTQ reads for base quality
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 22 -l 100 -o {output:q} <{input:q}"

# Kallisto quantification of one sample.
# Modified to declare the whole directory as the output.
rule kallisto_quant:
    output: directory("kallisto.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
       r"""mkdir {output:q}
           kallisto quant -i {input.index} -o {output:q} {input.fq1:q} {input.fq2:q} >& {output:q}/kallisto_quant.log
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
       r"""fastqc -o . {input:q}
           mv {wildcards.myfile:q}_fastqc.html {output.html:q}
           mv {wildcards.myfile:q}_fastqc.zip  {output.zip:q}
        """

rule salmon_quant:
    output: directory("salmon.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "salmon quant -i {input.index} -l A -1 {input.fq1:q} -2 {input.fq2:q} --validateMappings -o {output:q}"

rule salmon_index:
    output:
        idx = directory("{strain}.salmon_index")
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "salmon index -t {input.fasta} -i {output.idx} -k 31"

# A version of the MultiQC rule that ensures nothing unexpected is hoovered up by multiqc,
# by linking the files into a temporary directory.
rule multiqc:
    output:
        mqc_out = directory('multiqc_out'),
        mqc_in  = directory('multiqc_in'),
    input:
        salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
        fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
    shell:
       r"""mkdir {output.mqc_in:q}
           ln -snr -t {output.mqc_in:q} {input:q}
           multiqc {output.mqc_in:q} -o {output.mqc_out:q}
        """
