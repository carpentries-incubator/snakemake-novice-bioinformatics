echo "Test for episodes/files/ep07.Snakefile"

cp -vf episodes/files/ep07.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p multiqc

# Report should appear
test -s multiqc_out/multiqc_report.html

