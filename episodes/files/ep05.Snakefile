###
# Snakefile you should have after completing episode 05
###

# All conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# Alternative dynamic version using glob_wildcards, and with a print() to
# show us the content of the list:
#  CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition
#  print("Conditions are: ", CONDITIONS)

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
        log = "{strain}.kallisto_log",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {output.log}"
