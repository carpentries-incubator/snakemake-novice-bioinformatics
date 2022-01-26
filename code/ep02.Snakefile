###
# Snakefile you should have after completing episode 02
###

# Generic read counter rule using wildcards and placeholders
rule countreads:
  output: "{sample}.fq.count"
  input:  "reads/{sample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# Trim any FASTQ reads for base quality
rule trimreads:
  output: "trimmed/{sample}.fq"
  input:  "reads/{sample}.fq"
  shell:
    "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
