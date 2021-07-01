---
title: "Processing lists of inputs"
teaching: 0
exercises: 0
questions:
- "How do I process multiple files at once?"
- "How do I define a default set of outputs for my Snakefile?"
- "How do I make Snakemake decide what to process?"
- "How do I combine multiple files together?"
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
the command line. Often you want Snakemake to process all available samples. How can we do this?

The `yeast/reads` directory contains results from three conditions: `ref`, `etoh60` and `temp33`. There
are three replicates for each condition. There is minor inconsistency in the naming convention, as the `ref`
files do not have an underscore before the replicate number. Consistent naming is
important for Snakemake, so let's fix up the names before we go any further.

A good way to do this is by making symlinks, because that way you don't lose sight of the original file names.

~~~
$ mv reads original_reads
$ mkdir reads
$ ln -sr -t reads original_reads/*
$ for afile in reads/ref?_?.fq
> do
> newname=$( sed 's/ref/ref_/' <<<$afile )
> mv -v $afile $newname
> done
~~~
{: .language-bash}

These shell commands should symlink all the files and rename six of them. To tell Snakemake about our conditions
and replicates we can define some lists as Snakemake **global variables**.
Global variables can be added before the rules in the Snakefile.

~~~
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]
~~~

* The lists of quoted strings are enclosed in square brackets and comma-separated. If you know any Python you'll recognise
  this as Python list syntax.
* A good convention is to use capitalized names for these variables, but this is not mandatory.
* Although these are referred to as variables, you can't actually change the values once the workflow is running, so lists
  defined this way are more like constants.

## Using a Snakemake rule to define a batch of outputs

We'll add another rule to our Snakefile. This special target rule will have an `input` section but no `output` or `shell` sections (yet).

~~~
rule all_counts:
    input: expand("trimmed.{cond}_{rep}_1.fq.count", cond=CONDITIONS, rep=REPLICATES)
~~~

The `expand(...)` function in this rule generates a list of filenames, by taking the first thing in the parentheses as
a template and replacing `{cond}` with all the `CONDITIONS` and `{rep}` with all the `REPLICATES`. Since there are 3 of
each, this will yield 9 combinations - ie. 9 files we want to make.

This list goes into the `input` section of the rule. You might think that since these filenames are outputs that we want
from our workflow they should go into the `output` section. However, remember that `output`s of a rule are things the rule can
make itself, and this rule doesn't actually make anything. It's just a placeholder for a bunch of filenames.

We can now tell Snakemake to make all these files by using the target rule name on the command line, rather than the list of file
names:

~~~
$ snakemake -j1 -p all_counts
~~~

If you don't specify a target rule name or any file names on the command line when running Snakemake, the default is to use the
first rule in the Snakefile. So if `all_counts` is defined above the other rules you can simply say:

~~~
$ snakemake -j1 -p
~~~

> ## Adding a target rule to the Snakefile
>
> Check that the `all_counts` rule is working. Now adapt the Snakefile so that it makes all the counts for both of the pairs of
> reads (`_1.fq` and `_2.fq`), and also for both trimmed and untrimmed versions of the files. So you should end up with 36 count
> files in total.
>
>
> Tricky? Maybe a hint? Hmmm.
> > ## Solution
> >
> > ~~~
> > # Input conditions and replicates to process
> > CONDITIONS = ["ref", "etoh60", "temp33"]
> > REPLICATES = ["1", "2", "3"]
> > READ_ENDS  = ["1", "2"]
> > COUNT_DIR  = ["reads", "trimmed"]
> >
> > # Rule to make all counts at once
> > rule all_counts:
> >   input: expand("{indir}.{cond}_{rep}_{end}.fq.count", indir=COUNT_DIR, cond=CONDITIONS, rep=REPLICATES, end=READ_ENDS)
> > ~~~
> > {: .language}
> >
> > You can also put the lists directly into the `expand()` function if you think this is more readable. It's also possible to
> > split the function over more than one line, but note that you then need to start it on a new line too.
> >
> > ~~~
> > # Input conditions and replicates to process
> > CONDITIONS = ["ref", "etoh60", "temp33"]
> > REPLICATES = ["1", "2", "3"]
> >
> > # Rule to make all counts at once
> > rule all_counts:
> >   input:
> >     expand( "{indir}.{cond}_{rep}_{end}.fq.count", indir = ["reads", "trimmed"],
> >                                                    cond  = CONDITIONS,
> >                                                    rep   = REPLICATES,
> >                                                    end   = ["1", "2"] )
> > ~~~
> > {: .language}
> >
> {: .solution}
{: .challenge}

## Dynamically determining the inputs

At the start of this episode we ran a shell loop using a glob pattern: `reads/ref?_?.fq`. This matched all six files that we
needed to rename, and we didn't need to list them out explicitly. Snakemake allows you to do something similar with the
`glob_wildcards()` function.

~~~
CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition
print("Conditions are: ", CONDITIONS)
~~~

Here, the list of conditions is being set based on the files seen in the reads directory. The pattern used in glob_wildcards
looks much like the input and output parts of rules, with a wildcard in `{curly brackets}`, but here the pattern is being used
to search for matching files. We're only looking for read 1 of replicate 1 so this will return just three matches, and there is
just one wildcard in the glob pattern so we can assign this directly to a list. The `print()` statement will print out the value
of `CONDITIONS` when the Snakefile is run, and reassures us that the list really is the same as before.

~~~
$ snakemake -j1 -p all_counts
Conditions are:  ['etoh60', 'temp33', 'ref']
Building DAG of jobs...
Job counts:
	count	jobs
	1	all_counts
	36	countreads
	18	trimreads
	55
...
~~~

FIXME - I think this next bit can be converted into an exercise. Keep with the Python interpreter but ask some questions.

Using glob_wildcards() gets a little more tricky when you need a more complex match. To refine the match we can quickly test
out results by activating the Python interpreter.
This saves editing the Snakefile and running Snakemake just to see what `glob_wildcards()` will match. The Python interpreter is
like a special shell for Python commands, and because Snakemake functions are actually Python functions we can test them here.

~~~
$ python3
>>> from snakemake.io import glob_wildcards, expand
>>> glob_wildcards("reads/{condition}_1_1.fq")
Wildcards(condition=['etoh60', 'temp33', 'ref'])
~~~

This is the result we got before. So far, so good. Now to try listing the samples.

~~~
>>> glob_wildcards("reads/{sample}_1.fq")
Wildcards(sample=['temp33_3', 'temp33_2', 'etoh60_1', 'etoh60_3', 'ref_2', 'temp33_1', 'etoh60_2', 'ref_1', 'ref_3'])
~~~

Yeah this can be an exercise...

Before leaving the Python interpreter, we can also test the `expand()` function.

~~~
>>> SAMPLES = glob_wildcards("reads/{sample}_1.fq").sample
>>> expand("{sample}_{end}.fq.count", sample=SAMPLES, end=["1","2"])
~~~

> ## Glob with multiple wildcards
>
> If there are two wildcards in the glob pattern, dealing with the result becomes a little more tricky.
> Unless you're a Python programmer you probably don't want to start writing code like this, and for most
> cases in Snakemake there is no need to.
>
> However, for completeness, here is one way to re-combine two wildcards using `expand()` and `zip`.
>
> ~~~
> SAMPLES = expand( "{condition}_{samplenum}",
>                   zip,
>                   **glob_wildcards( "reads/{condition}_{samplenum}_1.fq" )._asdict() )
> ~~~
{: .callout}


## Rules that combine multiple inputs

Our `all_counts` rule is a rule which takes a list of input files. The length of that list is variable, and we've shown it can even be
calculated as part of the workflow by using `glob_wildcards()`. If we want to perform some combining operation on the list of files,
we can add `output` and `shell` sections to this rule.

In typical workflows, the final steps will combine all the results together into some big report. For our final workflow we'll be
doing this with *MultiQC*, but as a simple first example, let's just concatenate all the count files. In the shell this would be:

~~~
$ cat file1.count file2.count file3.count ... > all_counts.txt
~~~

If the input to a rule is a list of files, the `{input}` placeholder in the shell command will simply expand to that whole list,
separated by spaces. So to concetenate all the input files we just say:

~~~
shell:
  "cat {input} > {output}"
~~~

You can also combine named inputs and list inputs - any named input can be list of files rather than just a single file. Now,
using the `{input.name}` placeholder in the shell command will expand to that full list.


## Exercises

1. We already have an all_counts rule. Make it so that the rule concatenates all the count files into a single output file
   named 'all_counts_concatenated.txt'.
1. Adapt the rule so that trimmed and untrimmed reads are treated separately. That is, the rule will now have two named
   inputs rather than one single list. The output will still be a single concatenated file, with the untrimmed counts first and
   then the trimmed ones.

> ## Rules that make multiple outputs
>
> If we can have rules that combine lists of files, can we do the opposite and have a rule that produces a list of outputs?
> The answer is yes, but the situation is not completely symmetrical. Remember that Snakemake works out the full list of
> input and output files to every job *before* it runs a single job in the workflow. For a combining step, Snakemake will
> know how many samples/replicates are being combined from the start. For a splitting step, it may or may not be possible
> to predict the number of output files in advance. For cases where you really do need to handle a dynamic list of outputs,
> Snakemake has things called "dynamic rules" and "checkpoint rules". In practise these are very rarely needed, so we'll
> not be covering them here in the course.
>
{: .callout}


# Move to later...
~~~
$ multiqc -o multiqc_report .
~~~

We can define a new rule like so:

Things to note:

* The output of multiqc is a directory
* Something about supplying filenames to MultiQC
* ?? Do I need to add -f flag?

{% include links.md %}

