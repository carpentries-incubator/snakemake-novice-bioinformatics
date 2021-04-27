---
title: "Processing lists of inputs"
teaching: 0
exercises: 0
questions:
- "How do I process multiple files at once?"
- "How do I make Snakemake decide what to process?"
- "How do I combine multiple outputs together?"
objectives:
- "Use Snakemake to filter and count the lines in a FASTQ file"
keypoints:
- "Add key points"
---

We're introducing 'expand' and 'glob_wildcards' functions.
We'll also need to introduce the idea of lists and variables. Maybe do variables earlier?

--

## Defining a list of samples to process

So far, we've told Snakemake what output files to generate by giving the names of the desired files on
the command line. However, often you want Snakemake to process all available samples. How can we do this?

The `yeast/reads` directory contains results from three conditions: `ref`, `etoh60` and `temp33`. There
are three replicates for each condition. There is minor inconsistency in the naming convention, as the `ref`
files do not have an underscore before the replicate number, while the other two do. Consistent naming is
important for Snakemake, so let's fix up the names before we go any further.

A good way to do this is by making symlinks, because that way you don't lose sight of the original file names.

~~~
$ for afile in reads/ref?_?.fq
> do
> newname=$( sed 's/ref/ref_/' <<<$afile )
> ln -srv $afile $newname
> done
~~~

This shell loop should make you 6 new symlinks. To tell Snakemake about our conditions and replicates we can define some
lists as Snakemake **global variables**. Global variables go before the rules in the Snakefile.

~~~
CONDITIONS = ['ref', 'etoh60', 'temp33']
REPLICATES = ['1', '2', '3']
~~~

* The lists of quoted strings are enclosed in square brackets and comma-separated. If you know any Pythomn you'll recognise
  this as Python list syntax.
* A good convention is to use capitalized names for my variables, but this is not mandatory.
* Although these are referred to as variables, you can't actually change the values once the workflow is running, so lists
  defined this way are more like constants.

## Using a Snakemake rule to define a batch of outputs

We'll add another rule to our Snakefile...

~~~
rule all_counts:
    input: expand("trimmed.{cond}_{rep}_1.fq.count", cond=CONDITIONS, rep=REPLICATES)
~~~

The `expand(...)` function in this rule generates a list of filenames, by taking the first thing in the parentheses as
a template and replacing `{cond}` with all the `CONDITIONS` and `{rep}` with all the `REPLICATES`. Since there are 3 of
each, this will yield 9 combinations - ie. 9 output files.

This list goes into the `input` section of the rule. You might think that since these filenames are outputs that we want
to generate they should go into the `output` section. However, remember that `output`s of a rule are things the rule can
make, and this rule doesn't actually make anything. It's just a placeholder for a bunch of filenames.

We can now tell Snakemake to make all these files by using the rule name on the command line:

~~~
$ snakemake -j1 -p all_counts
~~~

Excercises:

Check that the `all_counts` rule is working. Now adapt the Snakefile so that it makes all the counts for both of the pairs of
reads, and also for both trimmed and untrimmed versions of the files. So you should end up with 36 count files in total.

Tricky? Maybe a hint? Hmmm.

## Dynamically determining the inputs

At the start of this episode we ran a shell loop using a glob pattern: `reads/ref?_?.fq`. This matched all six files that we
needed to symlink, and we didn't need to list them out explicitly. Snakemake allows you to do somethign similar with the
`glob_wildcards` function.

Do the glob here. Show that the result is the same.

{% include links.md %}

