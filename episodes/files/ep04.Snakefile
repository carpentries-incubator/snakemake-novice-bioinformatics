###
# Snakefile you should have after completing episode 04
# It's not modified during episode 05
###

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
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"

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
