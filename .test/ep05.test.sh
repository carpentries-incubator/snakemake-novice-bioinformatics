echo "Test for episodes/files/ep05.Snakefile"

cp -vf episodes/files/ep05.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# The all_counts rule assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p trimmed_counts_concatenated.txt kallisto.etoh60_1/abundance.h5

test -s trimmed_counts_concatenated.txt
test -s untrimmed_counts_concatenated.txt
test -s kallisto.etoh60_1/abundance.tsv
