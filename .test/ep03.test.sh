echo "Test for episodes/files/ep03.Snakefile"

cp -vf episodes/files/ep03.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

snakemake -j1 -p kallisto.etoh60_1/abundance.h5 trimmed.etoh60_2_2.fq.count

test -s trimmed/etoh60_2_2.fq
test -s kallisto.etoh60_1/abundance.h5
