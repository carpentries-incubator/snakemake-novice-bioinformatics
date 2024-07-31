echo "Test for episodes/files/ep06.Snakefile"

cp -vf episodes/files/ep06.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# This Snakefile assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p trimmed_counts_concatenated.txt \
                 trimmed.etoh60_1_1_fastqc.zip \
                 reads.etoh60_1_1_fastqc.zip

test -s trimmed_counts_concatenated.txt
test -s trimmed.etoh60_1_1_fastqc.html
test -s reads.etoh60_1_1_fastqc.html
