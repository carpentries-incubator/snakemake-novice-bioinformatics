###
# Snakefile you should have at the start of Episode 03
###

# Generic read counter rule using wildcards and placeholders
rule countreads:
  output: "{asample}.fq.count"
  input:  "reads/{asample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# Trim any FASTQ reads for base quality
rule trimreads:
  output: "trimmed/{asample}.fq"
  input:  "reads/{asample}.fq"
  shell:
    "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
