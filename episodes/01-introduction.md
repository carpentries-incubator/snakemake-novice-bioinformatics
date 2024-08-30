---
title: Running commands with Snakemake
teaching: 30
exercises: 30
---

::::::::::::::::::::::::::::::::::::::: objectives

- Create a Snakemake recipe (a Snakefile)
- Use Snakemake to count the lines in a FASTQ file

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I run a simple command with Snakemake?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Looking at the sample data

You should have the sample data files unpacked already (if not, refer back to the lesson setup).
Under `yeast/reads` you'll see a number of FASTQ files (extension `.fq`). These represent short
paired RNA-seq reads from an experiment on yeast cultured under three experimental conditions.
For now we'll just look at one single file, `ref1_1.fq`.

In the terminal:

```bash
$ cd snakemake_data/yeast
$ ls reads

$ head -n8 reads/ref1_1.fq
```

Files in FASTQ format have 4 lines for each sequence.

- Header line, beginning with `@`
- Sequence
- `+`
- Encoded quality per base

Let's count the number of lines in the file with `wc -l`. Again, in the shell:

```bash
$ wc -l reads/ref1_1.fq
58708
```

We can redirect this result to a file like so:

```bash
$ wc -l reads/ref1_1.fq > ref1_1.fq.count
$ head -v *.count
==> ref1_1.fq.count <==
58708 reads/ref1_1.fq
```

:::::::::::::::::::::::::::::::::::::::::  callout

## The sample dataset

The sample dataset represents a transcriptomics experiment in brewer's yeast
(Saccharomyces cerevisiae) under three conditions:

- **etoh60** - Treated with 60% ethanol
- **temp33** - Treated at 33 degrees celsius
- **ref**    - Untreated

For each condition there are 3 repeats, making 9 total samples. For each, total RNA
(or rather, cDNA) was sequenced on an Illumina HiSeq instrument.
For each repeat there is a pair of files as the sequencing is double-ended, so for example
`reads/etoh60_3_2.fq` contains the second of the read pairs from the third
ethanol-treated sample.

Don't worry about the biological meaning of this set-up. In the course, we're only going to get
as far as assessing the quality of the data and a very preliminary analysis. The data has been
subsampled to a fraction of the original size to make filtering and alignment operations fast,
so any deeper analysis is not going to yield meaningful results.

As well as the reads, we have a transcriptome for the yeast. This comes in a single FASTA file
(as opposed to the cDNA reads which are in FASTQ format) which has been GZip compressed.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Making a Snakefile

Within the `yeast` directory, edit a new text file named `Snakefile`.

Contents of `Snakefile`:

```source
rule countlines:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "wc -l reads/ref1_1.fq > ref1_1.fq.count"
```

:::::::::::::::::::::::::::::::::::::::  checklist

## Key points about this file

1. The file is named `Snakefile` - with a capital `S` and no file extension.
2. Some lines are indented. Indents must be with space characters, not tabs. See the
  [setup section](../learners/setup.md) for how to make your text editor do this.
3. The rule definition starts with the keyword `rule` followed by the rule name, then a colon.
4. We named the rule `countreads`. You may use letters, numbers or underscores, but the rule name
  must begin with a letter and may not be a keyword.
5. The keywords `input`, `output`, `shell` are all followed by a colon.
6. The file names and the shell command are all in `"quotes"`.
7. The output filename is given before the input filename. In fact, Snakemake doesn't care what
  order they appear in but we give the output first throughout this course. We'll see why soon.

::::::::::::::::::::::::::::::::::::::::::::::::::

Back in the shell we'll run our new rule. At this point, if there were any missing quotes, bad
indents, etc. we may see an error.

```bash
$ snakemake -j1 -F -p ref1_1.fq.count
```

For these early examples, we'll always run Snakemake with the `-j1`, `-F` and `-p` options. Later
we'll look more deeply at these and other available command-line options to Snakemake.

::::::::::::::::::::::::::::::::::: instructor

## Use of the -F flag

In the first few episodes we always run Snakemake with the `-F` flag, and it's not explained what
this does until Ep. 04. The rationale is that the default Snakemake behaviour when pruning the DAG
leads to learners seeing different output (typically the message "nothing to be done") when
repeating the exact same command. This can seem strange to learners who are used to scripting and
imperative programming.

The internal rules used by Snakemake to determine which jobs in the DAG are to be run, and which
skipped, are pretty complex, but the behaviour seen under `-F` is much more simple and consistent;
Snakemake simply runs every job in the DAG every time. You can think of `-F` as disabling the lazy
evaluation feature of Snakemake, until we are ready to properly introduce and understand it.

:::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Running Snakemake

Run `snakemake --help | less` to see the help for all available options.
What does the `-p` option in the `snakemake` command above do?

1. Protects existing output files
2. Prints the shell commands that are being run to the terminal
3. Tells Snakemake to only run one process at a time
4. Prompts the user for the correct input file

*Hint: you can search in the text by pressing `/`, and quit back to the shell with `q`*

:::::::::::::::  solution

## Solution

(2) Prints the shell commands that are being run to the terminal

This is such a useful thing we don't know why it isn't the default! The `-j1` option is what
tells Snakemake to only run one process at a time, and we'll stick with this for now as it
makes things simpler. The `-F` option tells Snakemake to always recreate output files, and
we'll learn about protected outputs much later in the course. Answer 4 is a total red-herring,
as Snakemake never prompts interactively for user input.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Counting sequences in a FASTQ file

FASTQ files contain 4 lines per sequence, as noted above. We can get BASH to divide the output of
`wc` by 4 to get the number of sequences:

```output
$ echo $(( $(wc -l <reads/ref1_1.fq) / 4 ))
14677
```

Note that the input filename is now preceeded by `<`. This is a little trick to get `wc` to print
only the number of lines, and not the filename. The `$( ... )` syntax captures this output value
and the `$(( ... ))` syntax encloses an arithmetic expression, which needs to be printed with
`echo`. Don't worry if this is unfamiliar - you just need to know that this is a shell command you
can copy and use to count the sequences.

:::::::::::::::::::::::::::::::::::: instructor

## Use of redirection (<) and shell arithmetic

A command is presented to count the sequences in a FASTQ file:

```
$ echo $(( $(wc -l <file.fq) / 4 ))
```

Understanding this in depth involves some advanced shell concepts that learners will not
necessarily be familiar with. However, other alternatives involve extra software, or use curly
brackets (which would have to be doubled-up) or are not robust (eg. `grep '^@' file.fq`).

:::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Counting sequences in Snakemake

Modify the Snakefile to count the number of **sequences** in `reads/ref1_1.fq`, rather than the
number of **lines**.

- Rename the rule to "countreads"
- Keep the output file name the same
- Remember that the result needs to go into the output file, not just be printed on the screen.

:::::::::::::::  solution

## Solution 1

```source
rule countreads:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "echo $(( $(wc -l <reads/ref1_1.fq) / 4 )) > ref1_1.fq.count"
```

:::::::::::::::::::::::::

Add a second rule to count the sequences in `reads/etoh60_1_1.fq`. Add this to the same Snakefile
you already made, under the "countreads" rule, and run your rules in the terminal. When running
the `snakemake` command you'll need to tell Snakemake to make both the output files.

:::::::::::::::  solution

## Solution 2

You can choose whatever name you like for this second rule, but it can't be "countreads" as
rule names need to be unique within a Snakefile. So in this example answer we
use "countreads2".

```source
rule countreads2:
  output: "etoh60_1_1.fq.count"
  input:  "reads/etoh60_1_1.fq"
  shell:
    "echo $(( $(wc -l <reads/etoh60_1_1.fq) / 4 )) > etoh60_1_1.fq.count"
```

Then in the shell...

```bash
$ snakemake -j1 -F -p ref1_1.fq.count etoh60_1_1.fq.count
```

If you think writing a separate rule for each output file is silly, you are correct. We'll
address this next!

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the full Snakefile](files/ep01.Snakefile) we should have by the end of
this episode.*



:::::::::::::::::::::::::::::::::::::::: keypoints

- Before running Snakemake you need to write a Snakefile
- A Snakefile is a text file which defines a list of rules
- Rules have inputs, outputs, and shell commands to be run
- You tell Snakemake what file to make and it will run the shell command defined in the
  appropriate rule

::::::::::::::::::::::::::::::::::::::::::::::::::


