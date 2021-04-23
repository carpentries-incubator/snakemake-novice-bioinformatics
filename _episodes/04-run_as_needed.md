---
title: "How snakemake runs only what is needed"
teaching: 0
exercises: 0
questions:
- "How does Snakemake avoid unecessary work?"
- "How do I control what steps will be run?"
objectives:
- "Understand the logic Snakemake uses when running and re-running rules"
keypoints:
- "Snakemake is awesome"
---

First, can I run Kallisto on the filtered reads because that makes the demo here easier? If not
I can introduce my re-pair script tho. OK, yes it works. Or at least it doesn't give an error and
the output looks like output.

## Snakemake is lazy and laziness is good

For the last few episodes, we've told you to run Snakemake like this:

~~~
$ snakemake -j1 -F -p desired_output_file
~~~
{: .language-bash}

As a reminder, the `-j1` flag tells Snakemake to run one job at a time, and `-p` is to print out
the shell commands before running them. Here we'll look at the `-F` flag which turns on `forceall`
mode.

At the end of the last chapter, we generated some kallisto results by running:

~~~
$ snakemake -j1 -F -p kallisto.ref1/abundance.h5
~~~
{: .language-bash}

Now try without the `-F` option. Assuming that the output files are already created, you'll see this:

~~~
$ snakemake -j1 -p kallisto.ref1/abundance.h5
Building DAG of jobs...
Nothing to be done.
Complete log: /home/zenmaster/data/yeast/.snakemake/log/2021-04-23T172441.519500.snakemake.log
~~~
{: .language-bash}

In normal operation, Snakemake only runs a job if:

1. An output file you explicitly requested is missing
1. The output of the job is missing and it is needed for another job
1. There is an input file which is newer than the output file

Let's demonstrate each of these in turn, by altering some files and re-running Snakemake.

rm kallisto.ref1/*

Now just re-runs the alignment

rm trimmed/ref1_?.fq

Nothing to be done - the output is missing but Snakemake already has the file you are telling it to
make so it doesn't worry.

touch transcriptome/xxx.gz

Describe what touch does.

Now Snakemake sees that one of the input files used in the whole process is newer than kallisto.ref1/abundance.h5,
so it will run the `kallisto index` and `kallisto quant` steps again. Of course, the `kallisto quant` step needs the
trimmed reads, so now the trimming step is re-run.

This behaviour is really useful when you want to:

1. Add new inputs to an existing analysis without re-processing everything
1. Continue running a workflow that failed part-way

In this case, the default Snakemake behaviour will just do the right thing. There are a few gotchas:

### Changing the rules

What if the rules change, rather than the input. Like, if you changed the filtering cutoffs?

### Removing input files

Snakemake can detect if you have added new input files, but not if you have removed input files. We'll look into this
more when we write rules with lists of files as input

### Incomplete jobs

Snakemake has a feature that it keeps a log of running jobs. If Snakemake crashes or exits unleanly then the next time
it runs it will refuse to use output files from incomplete jobs as the files may be partial output. We'll not look into
this during the course, but just be aware that this is what the Snakemake docs mean when they talk about incomplete jobs.
You tend to come across these more when working on a compute cluster, as opposed to a single machine.
{. :callout}

TODO - now convert stuff from the existing slides. And see what we can do regarding the exercises.
