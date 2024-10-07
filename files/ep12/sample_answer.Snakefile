# Run this Snakefile with the --use-conda flag to ensure the correct applications are present.

CONDITIONS = ["ref", "etoh60", "temp33"]
KMERS      = ["19", "21", "25", "27"]

rule all_max:
    input:
        expand("assem/{condition}_k{kmer}_max_contig.txt", condition = CONDITIONS,
                                                           kmer = KMERS )

rule cutadapt:
    output:
        read1 = "cutadapt/{sample}_1.fq",
        read2 = "cutadapt/{sample}_2.fq"
    input:
        read1 = "reads/{sample}_1.fq",
        read2 = "reads/{sample}_2.fq"
    params:
        adapter = "AGATCGGAAGAGC"
    conda: "assembly_conda_env.yaml"
    shell:
        r"""cutadapt -a {params.adapter} -A {params.adapter} \
            -o {output.read1} -p {output.read2} \
               {input.read1}     {input.read2}
         """

rule concatenate:
    output:
        read1 = "concatenate/{condition}_1.fq",
        read2 = "concatenate/{condition}_2.fq"
    input:
        read1s = ["cutadapt/{condition}_1_1.fq", "cutadapt/{condition}_2_1.fq", "cutadapt/{condition}_3_1.fq"],
        read2s = ["cutadapt/{condition}_1_2.fq", "cutadapt/{condition}_2_2.fq", "cutadapt/{condition}_3_2.fq"],
    shell:
        r"""cat {input.read1s} > {output.read1}
            cat {input.read2s} > {output.read2}
         """

rule assemble:
    output: "assem/{sample}_k{kmer}_contigs.fa"
    input:
        read1 = "concatenate/{sample}_1.fq",
        read2 = "concatenate/{sample}_2.fq"
    params:
        tmpdir = "velvet_tmp_{sample}_{kmer}"
    conda: "assembly_conda_env.yaml"
    shell:
        r"""velveth {params.tmpdir} {wildcards.kmer} -shortPaired -fastq -separate {input.read1} {input.read2}
            velvetg {params.tmpdir}
            mv {params.tmpdir}/contigs.fa {output}
        """

rule max_contig:
    output: "assem/{assem}_max_contig.txt"
    input:  "assem/{assem}_contigs.fa"
    conda: "assembly_conda_env.yaml"
    shell:
        "stats.sh {input} | grep 'Max contig length:' > {output}"

