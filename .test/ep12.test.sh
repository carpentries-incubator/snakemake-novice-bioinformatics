echo "Test for episodes/files/ep12.Snakefile"

cp -vf episodes/files/ep12.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

# Use dry-run to save a little time
res1=$( snakemake -Fn -p multiqc | \
          grep '^Would remove temporary output trimmed/.\+\.fq' | wc -l )
res2=$( snakemake -Fn -p multiqc | \
          grep '^Would remove temporary output .\+\.html' | wc -l )

# 9 samples so 18 trims
[[ "$res1" == 18 ]] && [[ "$res2" == 18 ]]

# Protected files don't show in the dry-run, so we actually have to make the Kallisto index
# and then check the perms on it. Then delete it so we don't upset later tests.
rm -f Saccharomyces_cerevisiae.R64-1-1.kallisto_index
snakemake -j1 -p Saccharomyces_cerevisiae.R64-1-1.kallisto_index
perms=$(stat --format %A Saccharomyces_cerevisiae.R64-1-1.kallisto_index)
[[ "$perms" == "-r--r--r--" ]]
rm -vf Saccharomyces_cerevisiae.R64-1-1.kallisto_index
