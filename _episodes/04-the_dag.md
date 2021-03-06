---
title: "How Snakemake plans what jobs to run"
teaching: 30
exercises: 15
questions:
- "How do I visualise a Snakemake workflow?"
- "How does Snakemake avoid unecessary work?"
- "How do I control what steps will be run?"
objectives:
- "View the DAG for our pipeline"
- "Understand the logic Snakemake uses when running and re-running rules"
keypoints:
- "Snakemake is awesome"
---

## The DAG

You may have noticed that one of the messages Snakemake always prints is:

~~~
Building DAG of jobs...
~~~

A DAG is a **Directed Acyclic Graph** and it can be pictured like so:

![DAG for our workflow][fig-dag]


The above DAG is based on our existing rules, and shows all the jobs Snakemake would run to trim, count and
quantify the "ref1" sample.

> ## Note that:
>
> * A rule can appear more than once, with different wildcards (a **rule** plus **wildcard values** defines a **job**)
> * The arrows show data flow, as well as dependency ordering between jobs
> * Snakemake can run the jobs in any order that doesn't break dependency - for example `kallisto quant` cannot run until
>   both `kallisto index` and `trimreads` have completed, but it may run before or after `countreads`
> * This is a work list, *not a flowchart*. There are no if/else decisions or loops - Snakemake runs every job in the DAG
>   exactly once
> * The DAG depends not only on the Snakefile, but on the requested target outputs and the files already present
{: .checklist}

> ## Question
>
> If we asked Snakemake to run 'kallisto quant' on all three of the reference samples (ref1, ref2, ref3), how many jobs would
> that be in total?
>
> > ## Answer
> >
> > 10 in total: 3 * kallisto_quant // 6 * trim // 1 * kallisto_index // 0 * countreads
> >
> {: .solution}
{: .challenge}

## Snakemake is lazy and laziness is good

For the last few episodes, we've told you to run Snakemake like this:

~~~
$ snakemake -j1 -F -p desired_output_file
~~~
{: .language-bash}

As a reminder, the `-j1` flag tells Snakemake to run one job at a time, and `-p` is to print out
the shell commands before running them.

The `-F` flag turns on `forceall` mode, and in normal usage you probably don't want this.

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

1. A target file you explicitly requested to make is missing
1. An output file is missing and it is needed in the process of making a target file
1. Snakemake can see an input file which is newer than an output file

Let's demonstrate each of these in turn, by altering some files and re-running Snakemake without the
`-F` option.

~~~
$ rm kallisto.ref1/*
$ snakemake -j1 -p kallisto.ref1/abundance.h5
~~~
{: .language-bash}

This just re-runs the kallisto quantification - the final step.

~~~
$ rm trimmed/ref1_?.fq
$ snakemake -j1 -p kallisto.ref1/abundance.h5
~~~
{: .language-bash}

"Nothing to be done" - some intermediate output is missing but Snakemake already has the file you are telling it to
make, so it doesn't worry.

~~~
$ touch transcriptome/*.fa.gz
$ snakemake -j1 -p kallisto.ref1/abundance.h5
~~~

The `touch` command is a standard Linux command which sets the timestamp of the file, so now the transcriptome looks
to Snakemake as if it was just modified.

Snakemake sees that one of the input files used in the process of producing `kallisto.ref1/abundance.h5` is newer than
the existing output file, so it needs to run the `kallisto index` and `kallisto quant` steps again. Of course, the
`kallisto quant` step needs the trimmed reads which we deleted earlier, so now the trimming step is re-run also.

This behaviour is really useful when you want to:

1. Change or add some inputs to an existing analysis without re-processing everything
1. Continue running a workflow that failed part-way

In these case, the default Snakemake behaviour will just do the right thing. In other situations you may need to
explicitly tell Snakemake to re-run some steps. For example...

## Changing the rules

If rules in the Snakefile change, rather than input files, Snakemake won't see that the results are out-of-date.
For example, if we changed the quality cutoffs within the trimreads rule then Snakemake would not automatically
re-run thos rules, because it only checks that the output file is newer than the input file.

The `-R` flag allows you to explicitly tell Snakemake that a rule has changed and that all outputs from that rule
need to be re-evaluated.

```
$ snakemake -j1 -Rtrimreads -p kallisto.ref1/abundance.h5
```

(Note - for some reason this only works if there is no space between `-R` and `trimreads`. This may be a bug in the
version of Snakemake we're using)

## Visualising the DAG

Snakemake can draw a picture of the DAG, if you run it like this:

~~~
snakemake -f --dag kallisto.etoh60_1/abundance.h5 | display -visual StaticColor
~~~

![DAG for partial workflow][fig-dag2]

The boxes drawn with dotted lines indicate steps that are not to be run, as the output files are newer than the
input files.

> ## Challenge
>
> Run `kallisto_quant` on all three of the **etoh60** samples. Now edit the Snakefile so that the `qual_threshold`
> for `trimreads` is "22", rather than "20".
>
> How would you get Snakemake to update all three Kallisto results:
>
>   1) By using the -R flag
>
>   2) By using the -f flag
>
>   3) By using the touch command
>
>   4) By deleting input files
>
> Use the `--dag` option as shown above to check your answers.
>
> > ## Solution
> >
> > To make all the Kallisto results in the first place:
> >
> > ~~~
> > snakemake -j1 -p kallisto.etoh60_{1,2,3}/abundance.h5
> > ~~~
> >
> > The {1,2,3} syntax is expanded by the shell into the 3 file names. You could also type all three names in full.
> >
> > 1) `$ snakemake -Rtrimreads --dag kallisto.etoh60_{1,2,3}/abundance.h5 | dot | display"`
> >
> > 2) `$ snakemake -j1 -p -f trimmed/etoh60_*.fq kallisto.etoh60_{1,2,3}/abundance.h5 --dag | dot | display`
> >
> > 3) `$ touch reads/etoh60_*`
> >
> > 4) `$ rm -r trimmed/etoh60_*.fq kallisto.etoh60_*`
> >
> > In general, getting Snakemake to re-run things by removing files is a bad idea. Using the `-R` flag or `touch` is
> > simpler and more reliable. If in doubt, and if it will not be too time consuming, use `-F` to run the whole
> > workflow from scratch.
> >
> {: .solution}
>
{: .challenge}

[fig-dag]: ../fig/dag_1.svg
[fig-dag2]: ../fig/dag_2.png

{% include links.md %}
