###
# Snakefile you should have after completing episode 01
###

# Count reads by counting lines in the file then dividing by 4 (using shell syntax)
rule countreads:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "echo $(( $(wc -l <reads/ref1_1.fq) / 4 )) > ref1_1.fq.count"

rule countreads2:
    output: "etoh60_1_1.fq.count"
    input:  "reads/etoh60_1_1.fq"
    shell:
        "echo $(( $(wc -l <reads/etoh60_1_1.fq) / 4 )) > etoh60_1_1.fq.count"
