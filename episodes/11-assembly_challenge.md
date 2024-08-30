---
title: Constructing a whole new workflow
teaching: 30
exercises: 90
---

::::::::::::::::::::::::::::::::::::::: objectives

- Apply what you have learned so far in a new example
- Understand how to diagnose errors in a workflow

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I approach making a new workflow from scratch?
- How do I understand and debug the errors I see?

::::::::::::::::::::::::::::::::::::::::::::::::::

## A script to assemble the sequences

This episode comprises a longer challenge to build a workflow from scratch, using the Snakemake
concepts we've learned through the course. A partial solution is provided as a Bash script, which
we will re-write in Snakemake and then extend it to get a full solution.

The script performs a de-novo assembly of the RNA sequences in the demo dataset. Note that this
does not produce a biologically meaningful assembly, but nevertheless this is the basis for
a nice example workflow, so we will go ahead and try.

The following shell script [can also be downloaded here](files/ep11/assembly_script.sh):

```bash
#!/bin/bash

# This script needs the following Bioconda packages:
#   cutadapt
#   velvet
#   bbmap

#  Use cutadapt to remove adapters. Note that cutadapt works on the paired sequences
mkdir cutadapt
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_1_1.fq -p cutadapt/ref_1_2.fq reads/ref_1_1.fq reads/ref_1_2.fq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_2_1.fq -p cutadapt/ref_2_2.fq reads/ref_2_1.fq reads/ref_2_2.fq
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o cutadapt/ref_3_1.fq -p cutadapt/ref_3_2.fq reads/ref_3_1.fq reads/ref_3_2.fq

# Combine all the samples. We still want to keep the read pairs separate.
mkdir concatenate
cat cutadapt/ref_1_1.fq cutadapt/ref_2_1.fq cutadapt/ref_3_1.fq > concatenate/ref_1.fq
cat cutadapt/ref_1_2.fq cutadapt/ref_2_2.fq cutadapt/ref_3_2.fq > concatenate/ref_2.fq

# Assemble the concatenated files with Velvet. This is a 2-step process: velveth then velvetg
kmer_len=21
velveth velvet_tmp_ref_${kmer_len} ${kmer_len} -shortPaired -fastq -separate concatenate/ref_1.fq concatenate/ref_2.fq
velvetg velvet_tmp_ref_${kmer_len}
mv velvet_tmp_ref_${kmer_len}/contigs.fa contigs.fa

# Find the longest contig using the BBMap stats script
stats.sh contigs.fa | grep 'Max contig length:' > max_contig.txt
head -v max_contig.txt
```

## Steps in the assembly process

The following figure presents the steps involved in the above script, when run on the `ref`
samples. You will see it already looks a little like a Snakemake DAG.

![][fig-flow]

* *Cutadapt* removes adapter sequence from the reads. Normally this should not be present but
  any that is left in will hinder the assembly. The `AGATCGGAAGAGC` sequence relates to the library
  prep kit used when preparing this cDNA for sequencing.
* *Velvet* is a simple (and old) short read assembler which attempts to piece together the raw
  short reads into longer sequences, known as contigs. More powerful de-novo assemblers exist, but
  for the purpose of this lesson Velvet does what we need. Velvet handles the type of paired short
  reads provided in the sample data.
* `velveth` and `velvetg` are components of the Velvet software which must be run one after the
  other to obtain a file of contigs (*contigs.fa*).
* When running `velveth` one must specify a *kmer length* parameter. We're not going to go into the
  theory of what this does, but the final Snakemake pipeline will help us explore the effect of
  changing this value.
* `stats.sh` is a handy program provided by the *BBMap* package, used here to tell us the longest
  contig.
* Note that we run *Cutadapt* on each `ref` sample separately, but then make a *Velvet* assembly
  of all the `ref` reads in one go.

:::::::::::::::::::::::::::::::::::::::  challenge

## Running the Bash script

Run the script as-is. In order to make it work, you should make a new conda environment named
"assembly-env" with the three packages installed.

How long is the longest contig found in the assembly that is produced?

:::::::::::::::  solution

## Solution

The commands to build and activate the Conda env can be found in the "Conda Integration"
chapter.

```bash
$ conda create -n assembly-env --channel bioconda --channel conda-forge cutadapt velvet bbmap
$ conda activate assembly-env
```

Now, after saving the script above into a text file, make it executable and run it.

```bash
$ chmod +x assembly_script.sh
$ ./assembly_script.sh
$ cat max_contig.txt
Max contig length:        655
```

Don't worry if you get a slightly different number. We've not specified the exact version of
the Velvet assembler so the final assembly may be slightly different if the software has been
updated.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## General tips for designing a Snakemake workflow

### Design phase

- Start out with pen and paper. Sketch how many rules you will have and how you expect the DAG
  to look.
- List any parameters or settings that might need to be adjusted.
- Work out which rules (if any) aggregate or split inputs (these may need input functions).
- Get the `input` and `output` parts of your Snakemake rules right before worrying about the
  `shell` sections. Remember that Snakemake builds the DAG before actually running the `shell`
  commands so you should `--dryrun` the workflow before running it for real, or even before you
  put all the `shell` commands in.
- Snakemake starts matching rules from the final target file to the input files, so you might
  prefer to work on your rules in this same order, starting by writing the last rule first.
- Choose the names for your *input*s, *output*s and *{wildcards}* carefully to make your Snakefile
  as readable as possible.

### Debugging a Snakemake workflow

As with all programming jobs, you will inevitably see bugs and errors the first time you try to run
a new Snakefile, because it's very easy to make typos. Sometimes Snakemake will give you a very
useful and precise error report, but other times less so. When trying to diagnose a problem,
remember that you need to establish which phase of execution failed, then look at the most common
error causes for that phase.

- **Parsing phase failures *(before Snakemake says 'Building DAG of jobs')*:**
  - Syntax errors
    - Missing commas or colons
    - Unbalanced quotes or brackets
    - Wrong indentation
    - ...etc
  - Failure to evaluate expressions
    - Any Python logic or you added outside of rules
    - Problems in functions (eg. `expand()`) in rule inputs or outputs
  - Other problems with rule definitions
    - Invalid rule names or declarations
    - Invalid wildcard names
    - Mismatched wildcards
- **DAG building failures *(before Snakemake tries to run any jobs)*:**
  - Failure to determine target output files
  - Multiple rules can make the same output
  - No rule to make a target file
    - Sometimes it is not clear why Snakemake wants to make that particular file
  - Circular dependency found (violating the "acyclic" property of a DAG).
  - An output file is write-protected
- **DAG running failures *(--dry-run works as expected, but a real run fails)*:**
  - If a job fails, Snakemake will report an error, delete all output files for that job,
    and stop.
  - Steps can fail if:
    - A placeholder in {curlies} is not a valid variable.
    - Any shell commands return non-zero status.
    - You reference any *$shell\_variable* which has not been set.
    - Any expected output file is missing after the commands have run.

If you get syntax errors, look out for unexpected or missing text colouring in GEdit which might
hint at the problem. Punctuation errors are very common: brackets, quotes, colons, line indentation
and commas all need to be just right. If you can't see the problem, ask a demonstrator for help.
Often just having a second pair of eyes on your code is all that is needed to spot a problem.

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge - building a full Snakemake workflow

Design a Snakemake workflow based upon the above script which will assemble all reads from each
of the three conditions (ref, etoh60 and temp33) and will do so at four different kmer lengths
(19, 21, 25 and 27), so twelve assemblies in total. Use your workflow to discover the length of
the longest contig in each of the twelve assemblies.

:::::::::::::::  solution

## Solution

A sample solution to this exercise [is available here](files/ep11/sample_answer.Snakefile),
along with a suitable [Conda environment spec](files/ep11/assembly_conda_env.yaml). However, there
is no single "correct" answer, so don't worry if your approach looks different. Hopefully we
will have time in the course to look through and compare some different answers.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


[fig-flow]: fig/assembly_flow.svg {alt='A box-and-arrow representation of the cutadapt,
concatenation, and assembly steps in the above script, with the names of the six files from the
three ref samples in pairs at the top. Arrows come down from each pair of files into a
corresponding cutadapt box, then under this the arrows cross as they go from the three cutadapt
boxes to the two boxes labelled "concatenate". Under this is a single box labelled
"assembly (velvet)" and a final box at the bottom labelled "find longest contig".'}


:::::::::::::::::::::::::::::::::::::::: keypoints

- By now, you are ready to start using Snakemake for your own workflow tasks
- You may wish to replace an existing script with a Snakemake workflow
- Don't be disheartened by errors, which are normal; use a systematic approach to diagnose
  the problem

::::::::::::::::::::::::::::::::::::::::::::::::::


