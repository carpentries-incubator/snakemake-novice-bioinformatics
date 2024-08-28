echo "Test for episodes/files/shadow.Snakefile"

cp -vf episodes/files/shadow.Snakefile snakemake_data/Snakefile
cd snakemake_data

# And run it
snakemake -F -j1 -p {normal,minimal,full}_out.txt

ls {normal,minimal,full}_out.txt
