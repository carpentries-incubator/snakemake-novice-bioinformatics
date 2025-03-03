echo "Test for episodes/files/ep09.Snakefile"

cp -vf episodes/files/ep09.Snakefile snakemake_data/yeast/Snakefile
cp -vf episodes/files/ep09/config.yaml snakemake_data/yeast/
cd snakemake_data/yeast

# This assumes renames are applied
( cd reads ; rename -v -s ref ref_ ref?_?.fq )

# Extract the default config from the comment in the Snakefile itself
sed -n -e '/config.yaml contents/,/^$/s/^# //p' Snakefile > config_fromfile.yaml

# They should be the same!
diff -qbZ config.yaml config_fromfile.yaml

# Use dry-run to save a little time
res1=$( snakemake -Fn --config -p multiqc | grep ^trimreads | tail -n1 )

# 9 samples, 2 .fq files each
[[ "$res1" =~ trimreads.+18 ]]

# And with alternative --config, only 2 samples precessed
res2=$( snakemake -Fn --config conditions='["temp33"]' replicates='["2", "3"]' -p multiqc | \
        grep ^trimreads | tail -n1 )

[[ "$res2" =~ trimreads.+4 ]]

