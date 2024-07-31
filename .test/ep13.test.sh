echo "Test for episodes/files/ep13.Snakefile"

cp -vf episodes/files/ep13.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

# Just do a dry-run for now, to test for basic errors.
snakemake -Fn -p multiqc
