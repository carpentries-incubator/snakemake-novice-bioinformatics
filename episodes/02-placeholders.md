---
title: Placeholders and wildcards
teaching: 40
exercises: 30
---

::::::::::::::::::::::::::::::::::::::: objectives

- Use Snakemake to count the sequences in any file
- Understand the basic steps Snakemake goes through when running a workflow
- See how Snakemake deals with some errors

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I make a generic rule?
- How does Snakemake decide what rule to run?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the Snakefile](code/ep01.Snakefile) you should have to start
the episode.*

## Wildcards and placeholders

In the previous episode you wrote two rules to count the sequences in two files. These work, but
they are not a very efficient use of Snakemake. We have eighteen input files to process and we
don't want to write eighteen near-identical rules!

To make a more general-purpose rule we need **placeholders** and **wildcards**. Here is a new rule
that will count the sequences in **any** of the `.fq` files.

```source
# New generic read counter
rule countreads:
    output: "{myfile}.fq.count"
    input:  "reads/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Comments in Snakefiles

In the above code, the line beginning `#` is a comment line. Hopefully you are already in the
habit of adding comments to your own scripts. Good comments make any script more readable, and
this is just as true with Snakefiles.

::::::::::::::::::::::::::::::::::::::::::::::::::

As a reminder, here's the non-generic version from the last episode:

```source
# Original version
rule countreads:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "echo $(( $(wc -l <reads/ref1_1.fq) / 4 )) > ref1_1.fq.count"
```

The new rule has replaced explicit file names with things in `{curly brackets}`, specifically
`{myfile}`, `{input}` and `{output}`.

### `{myfile}` is a **wildcard**

Wildcards are used in the `input` and `output` lines of the rule to represent parts of filenames.
Much like the `*` pattern in the shell, the wildcard can stand in for any text in order to make up
the desired filename. As with naming your rules, you may choose any name you like for your
wildcards, so here we used `myfile`. If `myfile` is set to `ref1_1` then the new generic rule will
have the same inputs and outputs as the original rule. Using the same wildcards in the input and
output is what tells Snakemake how to match input files to output files.

If two rules use a wildcard with the same name then Snakemake will treat them as different entities

- rules in Snakemake are self-contained in this way.

### `{input}` and `{output}` are **placeholders**

Placeholders are used in the `shell` section of a rule, and Snakemake will
replace them with appropriate values - `{input}` with the full name of the input file, and
`{output}` with the full name of the output file -- before running the command.

If we had wanted to include the value of the `myfile` wildcard directly in the `shell` command we
could have used the placeholder `{wildcards.myfile}` but in most cases, as here, we just need the
`{input}` and `{output}` placeholders.

:::::::::::::::::::::::::::::::::::::::  challenge

## Running the general-purpose rule

Modify your Snakefile to incorporate the changes described above, using the wildcard and
input/output placeholders. You should resist the urge to copy-and-paste from this workbook, but
rather edit the file by hand, as this will stick better in your memory.

You should delete the now-redundant second rule, so your Snakefile should contain just one rule
named *countreads*.

Using this new rule, determine: how many reads are in the `temp33_1_1.fq` file?

:::::::::::::::  solution

## Solution

After editing the file, run the commands:

```bash
$ snakemake -j1 -F -p temp33_1_1.fq.count

$ cat temp33_1_1.fq.count
20593
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Choosing the right wildcards

Our rule puts the sequence counts into output files named like `ref1_1.fq.count`. How would you
have to change the "countreads" rule definition if you wanted:

1) the output file for `reads/ref1_1.fq` to be `counts/ref1_1.txt`?

2) the output file for `reads/ref1_1.fq` to be `ref1_counts/fq.1.count` (for `reads/ref1_2.fq`
  to be `ref1_counts/fq.2.count`, etc.)?

3) the output file for `reads/ref1_1.fq` to be `countreads_1.txt`?

:::::::::::::::  solution

> ## Solution
> 
> In all cases, there is no need to change the `shell` part of the rule at all.
> 
> 1) > 
> ```source
> output: "counts/{myfile}.txt"
> input:  "reads/{myfile}.fq"
> ```
> 
> This can be done just by changing the `output:` line. You may also have considered the need to
> `mkdir counts` but in fact this is not necessary as Snakemake will create the output directory
> path for you before it runs the rule.
> 
> 2) > 
> ```source
> output: "{sample}_counts/fq.{readnum}.count"
> input:  "reads/{sample}_{readnum}.fq"
> ```
> 
> In this case, it was necessary to introduce a second wildcard, because the elements in the
> output file name are split up. The names chosen here are `{sample}` and `{readnum}` but you
> could choose any names as long as they match between the `input` and `output` parts. Once
> again, the output directory will be created for us by Snakemake, so the `shell` command does
> not need to change.
> 
> 3) This one **isn't possible**, because Snakemake cannot determine which input file you want to
>   count by matching wildcards on the file name "countreads\_1.txt". You could try a rule
>   like this:
> 
> ```source
> output: "countreads_{readnum}.count"
> input:  "reads/ref1_{readnum}.fq"
> ```
> 
> ...but it only works because "ref1" is hard-coded into the `input` line, and the rule will only
> work on this specific sample, not the other eight in our sample dataset. In general, input and
> output filenames need to be carefully chosen so that Snakemake can match everything up and
> determine the right input from the output filename.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Snakemake order of operations

We're only just getting started with some simple rules, but it's worth thinking about exactly what
Snakemake is doing when you run it. There are three distinct phases:

1. Prepares to run:
    1. Reads in all the rule definitions from the Snakefile
2. Plans what to do:
    1. Sees what file(s) you are asking it to make
    2. Looks for a matching rule by looking at the `output`s of all the rules it knows
    3. Fills in the wildcards to work out the `input` for this rule
    4. Checks that this input file is actually available
3. Runs the steps:
    1. Creates the directory for the output file, if needed
    2. Removes the old output file if it is already there
    3. Only then, runs the shell command with the placeholders replaced
    4. Checks that the command ran without errors *and* made the new output file as expected

For example, if we now ask Snakemake to generate a file named `wibble_1.fq.count`:

```output
$ snakemake -j1 -F -p wibble_1.fq.count
Building DAG of jobs...
MissingInputException in line 1 of /home/zenmaster/data/yeast/Snakefile:
Missing input files for rule countreads:
reads/wibble_1.fq
```

Snakemake sees that a file with a name like this could be produced by the *countreads* rule.
However, when it performs the wildcard substitution it sees that the input file would need to be
named `reads/wibble_1.fq`, and there is no such file available. Therefore Snakemake stops and gives
an error before any shell commands are run.

:::::::::::::::::::::::::::::::::::::::::  callout

## Dry-run (-n) mode

It's often useful to run just the first two phases, so that Snakemake will plan out the jobs to
run, and print them to the screen, but never actually run them. This is done with the `-n`
flag, eg:

```bash
$ snakemake -n -F -p temp33_1_1.fq.count
```

We'll make use of this later in the course.

::::::::::::::::::::::::::::::::::::::::::::::::::

The amount of checking may seem pedantic right now, but as the workflow gains more steps this will
become very useful to us indeed.

:::::::::::::::::::::::::::::::::::::::  challenge

## Adding a second rule

Here is a command that will trim and filter low quality reads from a FASTQ file.

```bash
$ fastq_quality_trimmer -t 20 -l 100 -o output.fq <input.fq
```

Add a second rule to your Snakefile to run this trimmer. You should make it so that valid outputs
are files with the same name as the input, but in a subdirectory named 'trimmed', for example:

- trimmed/ref1\_1.fq
- trimmed/temp33\_1\_1.fq
- *etc.*

:::::::::::::::  solution

## Solution

```source
# Trim any FASTQ reads for base quality
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
```

Bonus points if you added any comments to the code!

And of course you can run your new rule as before, to make one or more files at once. For
example:

```bash
$ snakemake -j1 -F -p trimmed/ref1_1.fq trimmed/ref1_2.fq
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## About fastq\_quality\_trimmer

`fastq_quality_trimmer` is part of the [FastX toolkit](https://hannonlab.cshl.edu/fastx_toolkit/)
and performs basic trimming on single FASTQ files.
The options `-t 20 -l 100` happen to be reasonable quality cutoffs for this dataset. This program
reads from standard input so we're using `<` to specify the input file, and the `-o` flag
specifies the output name.

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is a Snakefile](code/ep02.Snakefile) incorporating the changes made in
this episode.*



:::::::::::::::::::::::::::::::::::::::: keypoints

- Snakemake rules are made generic with placeholders and wildcards
- Snakemake chooses the appropriate rule by replacing wildcards such the the output matches the target
- Placeholders in the shell part of the rule are replaced with values based on the chosen wildcards
- Snakemake checks for various error conditions and will stop if it sees a problem

::::::::::::::::::::::::::::::::::::::::::::::::::


