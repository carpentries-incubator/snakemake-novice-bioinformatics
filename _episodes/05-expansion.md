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
the command line. Often you want Snakemake to process all available samples. How can we do this?

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
{: .language-bash}

This shell loop should make you 6 new symlinks. To tell Snakemake about our conditions and replicates we can define some
lists as Snakemake **global variables**. Global variables can be set before the rules in the Snakefile.

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

We'll add another rule to our Snakefile. This special target rule will have an `input` section but no `output` or `shell` sections.

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

> Excercises
>
> Check that the `all_counts` rule is working. Now adapt the Snakefile so that it makes all the counts for both of the pairs of
> reads (`_1.fq` and `_2.gq`, and also for both trimmed and untrimmed versions of the files. So you should end up with 36 count
> files in total.
>
>
> Tricky? Maybe a hint? Hmmm.
~~~
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]
READ_ENDS  = ["1", "2"]
COUNT_DIR  = ["reads", "trimmed"]

# Rule to make all counts at once
rule all_counts:
    input: expand("{indir}.{cond}_{rep}_{end}.fq.count", indir=COUNT_DIR, cond=CONDITIONS, rep=REPLICATES, end=READ_ENDS)
~~~

You can also put the lists directly into the `expand()` function if you think this is more readable. It's also possible to
split the function over more than one line, but note that you then need to start it on a new line too.

~~~
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# Rule to make all counts at once
rule all_counts:
    input:
        expand( "{indir}.{cond}_{rep}_{end}.fq.count", indir = ["reads", "trimmed"],
                                                       cond  = CONDITIONS,
                                                       rep   = REPLICATES,
                                                       end   = ["1", "2"] )
~~~

## Dynamically determining the inputs

At the start of this episode we ran a shell loop using a glob pattern: `reads/ref?_?.fq`. This matched all six files that we
needed to symlink, and we didn't need to list them out explicitly. Snakemake allows you to do something similar with the
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

Using glob_wildcards gets a little more tricky when you need a more complex match. For example, what if we want to list all the
samples rather than all the conditions? Firstly, we're going to quickly test out results by activating the Python interpreter.
This saves editing the Snakefile and running Snakemake just to see what `glob_wildcards` will match. The Python interpreter is
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
Wildcards(sample=['temp33_3', 'temp33_2', 'etoh60_1', 'etoh60_3', 'ref_2', 'temp33_1', 'ref3', 'etoh60_2', 'ref_1', 'ref_3', 'ref2', 'ref1'])

>>> glob_wildcards("reads/{condition}_{samplenum}_1.fq")
Wildcards(condition=['temp33', 'temp33', 'etoh60', 'etoh60', 'ref', 'temp33', 'etoh60', 'ref', 'ref'], samplenum=['3', '2', '1', '3', '2', '1', '2', '1', '3'])
~~~

Both of these are working, but neither is ideal.  The second, using two wildcards, yields two lists. It's possible to
re-combine these lists to construct the list of samples, but it's tricky.

The first is picking up both the `ref1` and `ref_1` files because the wildcard matches either. We could avoid this by
cleaning up the input directory, but we can also be more specific about what the wildcards can match.

~~~
>>> glob_wildcards("reads/{sample,.+_.}_1.fq")
Wildcards(sample=['temp33_3', 'temp33_2', 'etoh60_1', 'etoh60_3', 'ref_2', 'temp33_1', 'etoh60_2', 'ref_1', 'ref_3'])
~~~

This does the trick. We added a regex pattern to the wildcard to restrict what it can match. Note the syntax:

* The wildcard name is folloed by a comma, then the regex pattern
* The expression is `.+_.` meaning:
 * `.+` One or more of any character
 * `_` Followed by a literal underscore
 * `.` Followed by a single character

We could have been even more prescriptive, but this does the job. Now we can put this into the Snakefile:

~~~
SAMPLES = glob_wildcards("reads/{sample,.+_.}_1.fq")
~~~

> # Info box
>
> For completeness, here's one way to work with two wildcards. I'll present this here without explanation, but
> you'll see it's not pretty.
>
> ~~~
> SAMPLES = expand( "{condition}_{samplenum}",
>                   zip,
>                   **glob_wildcards( "reads/{condition}_{samplenum}_1.fq" )._asdict() )
> ~~~

## Directories as inputs and outputs

In the exercise below, we'll work with a program that produces a whole directory of files as output. We already saw this
for `kallisto quant` and in this case the directory contained three files, so we listed these as three outputs of the rule.

~~~
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    ...
~~~

There are two other valid approaches:

1. List just a subset of the output files as outputs of the rule. Snakemake does not care (in fact does not check) if the
   command produces other files too. So for the `kallisto_quant` rule we could have just said `output: "abundance.h5"` and
   the rule would work. The other outputs still get crrated but Snakemake does not condifer them when linking rules.

1. Tell snakemake that the output of the rule is a directory. This way, Snakemake will not consider the files inside the
   directory at all. Do this by adding `directory(...)` around the output path.

~~~
rules kallisto_quant:
    output: directory("kallisto.{sample}")
    ...
~~~

Note that you only have to do this for outputs. The input to a rule may be a directory without the need for any special
syntax.

Exercise needed here. glob_wildcards gets hairy when used with mutiple wildcards since you get back a named tuple of lists.

> ## Challenge
>
> Adapt the Snakefile to run 'kallisto quant' on all 9 samples, that is all three repeats of all three conditions. Rather than
> listing the inputs explicitly, use `glob_wildcards` to find them and make a 
>
> An alternative to kallisto for transcript quantification is `salmon`. The procedure is virtually identical, having an indexing
> step and a quantification step. Note that in real usage one is advised to add decoy sequences to the index but for the
> purposes of this tutorial we'll just keep things simple.
>
> Based upon the following commands:
>
> ~~~
> $ salmon index -t <transcriptome as fastq> -i <index name> -k 31
> $ salmon quant -i <index name> -l A -1 <fastq1> -2 <fastq2> --validateMappings -o <output path>
> ~~~
>
> Add rules to build a transcriptome index and run quantification using salmon. Keep the existing Kallisto rules so that the pipeline
> now supports both options. Amend the target rule to run all 18 quantifications, that is all 9 samples in Salmon and all 9 in
> Kallisto.
>


CONDITIONS = glob_wildcards("reads/{condition}_{read}_1.fq").condition

{% include links.md %}

