echo "Test for episodes/files/ep03.Snakefile"

cp -vf episodes/files/ep03.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

snakemake -j1 -F -p ref1_1.reads_removed.txt etoh60_1_1.reads_removed.txt

head -v *_removed.txt

[[ $(<ref1_1.reads_removed.txt) == 399 ]]
[[ $(<etoh60_1_1.reads_removed.txt) == 532 ]]
