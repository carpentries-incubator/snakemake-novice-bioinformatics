# A Snakefile to demonstrate how shadow rules work. Try:
#
# $ snakemake -F -j1 -s shadow.Snakefile {normal,minimal,full}_out.txt
# $ head -n 30 {normal,minimal,full}_out.txt
# $ ls ?_temp_file

rule normal_rule:
    output: "normal_out.txt"
    shell:
       r"""exec >{output}
           echo no shadow mode
           echo Current directory is: `pwd`
           touch 1_temp_file
           tree
        """

rule minimal_shadow_rule:
    output: "minimal_out.txt"
    shadow: "minimal"
    shell:
       r"""exec >{output}
           echo minimal shadow mode
           echo Current directory is: `pwd`
           touch 2_temp_file
           tree
        """

rule full_shadow_rule:
    output: "full_out.txt"
    shadow: "full"
    shell:
       r"""exec >{output}
           echo full shadow mode
           echo Current directory is: `pwd`
           touch 3_temp_file
           tree
        """
