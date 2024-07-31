# Snakefile tests

To ensure that all the sample answers are really run-able, this directory contains test scripts for
all of the files. They are designed to be run from the .github/workflows/code-test.yml action upon
commits to GitHub but you can run the scripts locally as long as you unpack the `snakemake_data`
directory into the root of the GIT checkout (ie. where the episodes/files directory resides).
