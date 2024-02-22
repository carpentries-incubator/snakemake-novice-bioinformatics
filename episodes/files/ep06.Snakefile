###
# Snakefile you should have after completing episode 06
###

# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# Rule to make all counts and compile the results in two files
rule all_counts:
    input:
        untrimmed = expand( "reads.{cond}_{rep}_{end}.fq.count",   cond  = CONDITIONS,
                                                                   rep   = REPLICATES,
                                                                   end   = ["1", "2"] ),
        trimmed   = expand( "trimmed.{cond}_{rep}_{end}.fq.count", cond  = CONDITIONS,
                                                                   rep   = REPLICATES,
                                                                   end   = ["1", "2"] ),
    output:
        untrimmed = "untrimmed_counts_concatenated.txt",
        trimmed   = "trimmed_counts_concatenated.txt",
    shell:
        "cat {input.untrimmed} > {output.untrimmed} ; cat {input.trimmed} > {output.trimmed}"

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
# Modified to declare the whole directory as the output.
rule kallisto_quant:
    output: directory("kallisto.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o {output} {input.fq1} {input.fq2}"

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
