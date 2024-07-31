echo "Test for episodes/files/ep10.Snakefile"

cp -vf episodes/files/ep10.Snakefile snakemake_data/yeast/Snakefile
cd snakemake_data/yeast

# Emit salmon-1.2.1.yaml (learners are advised to install these versions then freeze the env)
cat >salmon-1.2.1.yaml <<.
name: salmon-1.2.1
channels:
  - conda-forge
  - bioconda
dependencies:
  - salmon=1.2.1
  - tbb=2020.2
.

# Now run the thing
snakemake -j1 --sdm conda -p salmon.etoh60_1

# A conda env should have been created
conda env list | grep -F /.snakemake/

# But remove it
for e in $(conda env list | grep -F /.snakemake/) ; do echo conda env remove --yes -p "$e" ; done

# We can confirm which version of Salmon was run, according to salmon.etoh60_1/cmd_info.json
salmon_version=$( jq -r .salmon_version salmon.etoh60_1/cmd_info.json )
[[ "$salmon_version" = 1.2.1 ]]

