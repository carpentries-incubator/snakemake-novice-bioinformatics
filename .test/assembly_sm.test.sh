echo "Test for episodes/files/assembly_with_conda.Snakefile"

cp -vf episodes/files/ep11/sample_answer.Snakefile snakemake_data/yeast/Snakefile
cp -vf episodes/files/ep11/assembly_conda_env.yaml snakemake_data/yeast/

cd snakemake_data/yeast

# The Snakefile assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

# Use all the cores
snakemake -c all -p --sdm conda

# We should have made a new conda env, but remove it
for e in $(conda env list | grep -F /.snakemake/) ; do echo conda env remove --yes -p "$e" ; done

ls assem
# We should have 12 assemblies
[[ $(echo assem/*_max_contig.txt | wc -w) == 12 ]]
