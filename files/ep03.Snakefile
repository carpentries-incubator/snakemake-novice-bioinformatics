###
# Snakefile you should have after completing episode 03
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
