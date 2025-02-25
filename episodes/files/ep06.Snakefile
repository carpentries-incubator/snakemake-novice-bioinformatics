###
# Snakefile you should have after completing episode 06
###

# All conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# Alternative dynamic version using glob_wildcards, and with a print() to
# show us the content of the list:
#  CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition
#  print("Conditions are: ", CONDITIONS)

# Rule to make all counts, calculate the difference, and compile the results in two files
rule all_differences:
    input: expand("all_read{end}_removed.txt", end=["1","2"])

rule all_differences_per_end:
    input:
        expand( "{cond}_{rep}_{{end}}.reads_removed.txt", cond  = CONDITIONS,
                                                          rep   = REPLICATES )
    output:
        "all_read{end}_removed.txt"
    shell:
        "cat {input} > {output}"

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

# Find the difference between untrimmed and trimmed count files
rule calculate_difference:
    output: "{myfile}.reads_removed.txt"
    input:
        untrimmed = "reads.{myfile}.fq.count",
        trimmed = "trimmed.{myfile}.fq.count",
    shell:
        "echo $(( $(<{input.untrimmed}) - $(<{input.trimmed}) )) > {output}"

# Kallisto quantification of one sample
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"

rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz",
    log: "{strain}.kallisto_log"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {log}"
