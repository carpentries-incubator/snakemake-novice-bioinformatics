---
title: "Placeholders and Wildcards"
teaching: 0
exercises: 0
questions:
- "How do I make a generic rule?"
objectives:
- "Use Snakemake to count the sequences in any file"
- "Understand the basic steps Snakamake goes through when running a workflow"
- "See how Snakemake deals with missing inputs"
keypoints:
- "Add key points"
---

(Similar to the workflows-snakemake thing that already exists. It intros wildcards after getting the first
workflow running. I can see why this is a good idea.)

## Wildcards and placeholders

In the previous episode you wrote two rules to count the sequences in two files. These work, but they are not
a very efficient use of Snakemake. We have eighteen input files to process and we don't want to write eighteen
near-identical rules.

To make a more general-purpose rule we need **placeholders** and **wildcards**. Here is a new rule that will count the
sequences in **any** of the `.fq` files.

~~~
# New generic read counter
rule countreads:
  output: "{asample}.fq.count"
  input:  "reads/{asample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~
{: .language}

> ## Comments in Snakefiles
>
> In the above code, the line beginning `#` is a comment line. Hopefully you are already in the habit of adding
> comments to your own scripts. Good comments make any script more readable, and this is just as true with
> Snakefiles.
>
{: .callout}

As a reminder, here's the non-generic version from the last episode:

~~~
# Original version
rule countreads:
  output: "ref1_1.fq.count"
  input:  "reads/ref1_1.fq"
  shell:
    "echo $(( $(wc -l <reads/ref1_1.fq) / 4 )) > ref1_1.fq.count"
~~~
{: .language}

The new rule has replaced file names with things in `{curly brackets}`, specifically `{asample}`, `{input}` and `{output}`.

`{asample}` is a **wildcard**. Wildcards are used in the `input` and `output` lines of the rule to represent parts of filenames.
As with rule names, you may choose any name you like for your wildcards, so here we chose `asample`.
If two rules use a wildcard with the same name then Snakemake will treat them as completely different - rules in Snakemake
are self-contained in this way.

`{input}` and `{output}` are **placeholders**. Placeholders are used in the `shell` section of a rule, and Snakemake will
replace them with appropriate values - `{input}` with the full name of the input file, and `{output}` with the full name of
the output file -- before running the command.

If we had wanted to include the value of the `asample` wildcard directly in the `shell` command we sould have used the placeholder
`{wildcards.asample}` but in most cases, as here, we just need the `{input}` and `{output}` placeholders.

> ## Running the general-purpose rule
>
> Modify your Snakefile to incorporate the changes described above, using the wildcard and input/output placeholders.
> You should resist the urge to copy-and-paste from this workbook, but rather edit the file by hand, as this will stick
> better in your memory.
>
> You should delete the now-redundant second rule, so your Snakefile should contain just one rule named "countreads".
>
> Using this new rule, determine: how many reads are in the temp33_1_1.fq file?
>
> > ## Solution
> >
> > After editing the file, run the commands:
> >
> > ~~~
> > $ snakemake -j1 -F -p temp33_1_1.fq.count
> >
> > $ cat temp33_1_1.fq.count
> > #???#
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## Choosing the right wildcards
>
> Our rule puts the sequence counts into output files named like `ref1_1.fq.count`. How would you have to change the "countreads"
> rule definition if you wanted:
>
>  1) the output file for `reads/ref1_1.fq` to be `counts/ref1_1.txt`?
>
>  2) the output file for `reads/ref1_1.fq` to be `ref1_counts/fq.1.count` (for `reads/ref1_2.fq` to be `ref1_counts/fq.2.count`)?
>
>  3) the output file for `reads/ref1_1.fq` to be `countreads_1.txt`?
>
> > ## Solution
> >
> > 1)
> > ~~~
> > output: "counts/{asample}.txt"
> > input:  "reads/{asample}.fq"
> > ~~~
> > {: .language}
> >
> > This can be done just by changing the `output:` line. You may also have considered the need to `mkdir counts` but
> > in fact this is not necessary as Snakemake will create the output directory path for you before it runs the rule.
> >
> > 2)
> > ~~~
> > output: "{asample}_counts/fq.{readnum}.count"
> > input:  "reads/{asample}_{readnum}.fq"
> > ~~~
> > {: .language}
> >
> > In this case, it was necessary to introduce a second wildcard, because the elements in the output file name
> > are split up. The name chosen here is `{readnum}` but you could choose any name as long as the names match
> > between the `input` and `output` parts. Once again, the output directory will be created for us by Snakemake,
> > so the `shell` command does not need to change.
> >
> > 3) This one isn't possible, because Snakemake cannot determine which input file you want to count by matching wildcards
> > on the file name "countreads_1.txt". In general, input and output filenames need to be carefully chosen so that Snakemake
> > can match everything up. There's something of an art to this!
> >
> {: .solution}
>
{: .challenge}


## Snakemake order of operations

We're only just getting started with some simple rules, but it's worth thinking about exactly what Snakemake is doing
when you run it:

1. Prepares to run:
    1. Reads in all the rule definitions from the Snakefile
1. Plans what to do:
    1. Sees what file(s) you are asking it to make
    1. Looks for a matching rule by looking at the `output`s of all the rules it knows
    1. Fills in the wildcards to work out the `input` for this rule
    1. Checks that this input file is actually available
1. Runs the steps:
    1. Checks that the directory for the output file is available, and if not creates it
    1. Only then, runs the shell command with the placeholders replaced
    1. Checks that the command ran without errors *and* made the output file as expected

For example, if we now ask Snakemake to generate `wibble_1.fq.count`:

~~~
$ snakemake -j1 -F -p wibble_1.fq.count
Building DAG of jobs...
MissingInputException in line 1 of /home/zenmaster/data/yeast/Snakefile:
Missing input files for rule countreads:
reads/wibble_1.fq
~~~

Snakemake sees that a file with a name like this could be produced by the "countreads" rule. However, when it performs
the wildcard substitution it sees that the input file would need to be named `reads/wibble_1.fq`, and there is no such
file available. Therefore Snakemake stops and gives an error before running a single command.

The amount of checking may seem pedantic right now, but as the workflow gains more steps this will become very useful
to us indeed.

> ## Adding a second rule
>
> Here is a command that will trim and filter low quality reads from a FASTQ file.
>
> ~~~
> $ fastq_quality_trimmer -t 20 -l 100 -o output.fq <input.fq
> ~~~
> {: .language-bash}
>
> Add a second rule to your Snakefile to run this trimmer. You should make it so that valid outputs are files with
> the same name as the input, but in a subdirectory named 'trimmed', for example:
>
> * trimmed/ref1_1.fq
> * trimmed/temp33_1_1.fq
> * *etc.*
>
> > ## Solution
> >
> > ~~~
> > # Trim any FASTQ reads for base quality
> > rule trimreads:
> >   output: "trimmed/{asample}.fq"
> >   input:  "reads/{asample}.fq"
> >   shell:
> >     "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
> > ~~~
> > {: .language}
> >
> > Bonus points if you added any comments to the code!
> >
> > And of course you can run your new rule as before, to make one or more files at once. For example:
> >
> > ~~~
> > $ snakemake -j1 -F -p trimmed/ref1_1.fq trimmed/ref1_2.fq
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## About fastq_quality_trimmer
>
> `fastq_quality_trimmer` is part of the FastX toolkit and performs basic trimming on single FASTQ files.
> The options `-t 20 -l 100` happen to be reasonable quality cutoffs for this dataset. This program reads from
> standard input so we're using `<` to specify the input file, and the `-o` flag specifies the output name.
{: .callout}

More here??

Some multiple choice? Like what's wrong with this Snakefile?

Missing inputs - what if you ask for a file for which there is no rule? Fairly self-explanatory. Demo this.

What if you ask for a file and there is a rule, but no input? "No rule to make" error. Note that Snakemake
does not even try to run the rule.

What if the rule doesn't generate the output file? Give an example here. Or maybe ask it as a question.

Q: Break your rule so that it doesn't produce an output file (be more specific). What happens when you run it?

Hmmm. Maybe this not so good. We'll come back to it.

{% include links.md %}

