echo "Test for episodes/files/ep07.Snakefile"

cp -vf episodes/files/ep07.Snakefile snakemake_data/yeast/Snakefile
cp -vf episodes/files/ep07/wrappers.Snakefile snakemake_data/yeast/wrappers.Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p multiqc

# Report should appear
test -s multiqc_out/multiqc_report.html

# There is also the version with wrappers
snakemake -j1 --delete-all-output multiqc
test '!' -e multiqc_out

snakemake -j1 -s wrappers.Snakefile -p multiqc

# Report should appear again
test -s multiqc_out/multiqc_report.html
