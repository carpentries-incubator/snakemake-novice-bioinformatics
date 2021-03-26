---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "How do I run a simple command with Snakemake?"
objectives:
- "Use Snakemake to count the lines in a file"
keypoints:
- "Add key points"
---

Intro the .fq files. Note they are chicken RNA and remind the FASTQ format:

* Header line
* Sequence
* +
* Encoded quality

Let's count the number of lines in a file with 'wc -l'

  wc -l data/liver_1.fq

We can redirect this info to a file like so:

 wc -l data/liver_1.fq > liver_1.fq.count

Making a Snakefile:

rule: countlines
  input:  "data/liver_1.fq"
  output: "liver_1.fq.count"
  shell:
    "wc -l data/liver_1.fq > liver_1.fq.count"

Editor setup here with some explanation.

How to run the thing and some potential problems you might see (missing quotes, bad indents, missing input files, nothting to be done)

Exercises:

FASTQ files contain 4 lines per read. We can get BASH to do do the division for us to get the number of sequences:

echo $(( $(wc -l <somefile.fq) / 4 ))

Try this in the shell first. Modify the Snakefile to count the number of sequences in data/liver_1.fq.

Add another rule to count the sequences in gut_1.fq. If you this rule in the file above the countlines rule then Snakemake will run
it by default. If you put it after, you will need to explicitly tell snakemake the rule to be run.

{% include links.md %}

