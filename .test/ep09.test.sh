echo "Test for episodes/files/ep09.Snakefile"

cp -vf episodes/files/ep09.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

# Use dry-run to save a little time
res1=$( snakemake -Fn -p multiqc | grep 'kallisto quant -t 4' | wc -l )

# 9 samples so 9 calls to kallisto
[[ "$res1" == 9 ]]

# Check salmon as well
res2=$( snakemake -Fn -p multiqc | grep 'salmon quant -p 4' | wc -l )

[[ "$res2" == 9 ]]

