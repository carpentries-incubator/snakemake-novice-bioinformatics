echo "Test for episodes/files/ep02.Snakefile"

cp -vf episodes/files/ep02.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

snakemake -j1 -p trimmed/etoh60_1_1.fq etoh60_1_1.fq.count

test -s trimmed/etoh60_1_1.fq
test -s etoh60_1_1.fq.count
