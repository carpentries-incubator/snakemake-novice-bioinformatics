echo "Test for episodes/files/ep06.Snakefile"

cp -vf episodes/files/ep06.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# The all_counts rule assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p all_read1_removed.txt kallisto.etoh60_1/abundance.h5
test -s all_read1_removed.txt
test -s kallisto.etoh60_1/abundance.tsv
test '!' -s all_read2_removed.txt

snakemake -j1 -p all_differences
test -s all_read2_removed.txt
