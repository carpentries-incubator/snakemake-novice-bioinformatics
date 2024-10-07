echo "Test for episodes/files/ep08.Snakefile"

cp -vf episodes/files/ep08.Snakefile snakemake_data/yeast/Snakefile
cp -vf episodes/files/ep08/wrappers.Snakefile snakemake_data/yeast/wrappers.Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

snakemake -j1 -p multiqc

# Report should appear
test -s multiqc_out/multiqc_report.html

# There is also the version with wrappers
snakemake -j1 --delete-all-output multiqc
test '!' -e multiqc_out

echo "Testing the version that uses wrappers..."

snakemake -c all -s wrappers.Snakefile -p multiqc

# Report should appear again
test -s multiqc_out/multiqc_report.html
