---
title: 'Prerequisites'
---

## Prerequisites

This lesson assumes you have a working knowledge of running bioinformatics tools in the Bash shell.

Most of the core concepts are introduced in the *Unix Shell* carpentries core lesson, but the
hands-on experience of having analysed some of your own data within the Linux shell environment
is what will really prepare you for this course.

A detailed breakdown of the concepts and syntax that comes up in the course is given below.

## What you will not learn

If you are looking for guidance in how to select appropriate modern tools for your bioinformatics
analysis then you are in the wrong place. We deliberately use some older but simple and stable
tools like *fastx toolkit* and the *velvet* assembler. And we run our tools with simple defaults
rather than considering the details of proper configuration. The rationale is that we are
presenting these programs to illustrate how Snakemake works, not to learn about the tools and
approaches themselves.

Also, if you are working under pressure to analyse your data as soon as possible, your time may not
be best spent learning Snakemake. The ability to orchestrate analyses with a workflow system is an
essential ability for a modern bioinformatician, but if you are struggling to make your analysis
work with regular shell commands, scripts and web resources then adding Snakemake to the mix only
introduces another layer of complexity. Having said that, once you are confident using Snakemake
you will find it comes in handy even for your simple analysis tasks.

## Specific terms and concepts

### Editing code

Snakemake is a programming language and programming requires editing code in a text editor. There
are many many such editors available and most are suitable for use with Snakemake. The setup
instructions for this course have some pointers for the *GEdit* and *nano* editors, but the core
abilities you need are:

* Start the text editor
* Load and save a file
* Edit, copy, and paste text
* Turn on syntax highlighting in the editor
* Show line numbers in the editor

### Shell syntax and common commands

You should, seeing this example Bash shell command:

```
$ ls -ltrF --directory --human-readable foo_*
```

* See that `$` represents the *shell prompt*
* Know that `ls` is a command for listing files and directories
* Infer that `-ltrF` represent four single-letter options to the `ls` command
* Infer that `--directory` and `--human-readable` represent long-form opions
* Know how to look up the meaning of such options in the `man` page
* See that `foo_*` is a wildcard (aka. glob) pattern that may match multiple file names

And these specific concepts occur within the course:

* Use `cd` to change the current working directory
* Use `&` to run a command as a background job
* Use `>` and `<` to redirect terminal output to, and input from, files
* Use `|` to redirect (pipe) output directly between commands
* Use `wget` to fetch a file from a web URL
* Create and remove directories with `mkdir` and `rmdir`
* Use `rm` and `rm -r` to remove files and non-empty directories
* Use `ln -s` to create a symbolic link
* Use `cat` and `head` to show the contents of a text file in the terminal
* Run a Bash script file containing shell commands: `bash scriptname.sh`

## Biology and bioinformatics

This course is geared to bioinformaticians, and those studying DNA sequences, so those without
a background in this research may struggle to understand the motivation behind the example
workflows. There follows a very brief and simplified synopsis of the background concepts. You may
want to dig into these further, to which end some links to Wikipedia are provided.

[Second generation DNA sequencing](https://en.wikipedia.org/wiki/Massive_parallel_sequencing) is
a technology allowing large quantities of DNA to be read quickly on a single machine, by sequencing
millions of short DNA fragments in parallel. These fragments are typically read from
both ends, yielding a matched pair of sequence reads, but the technology can only see a short
way into the fragment from either end. If the physical fragments are short enough the pair of reads
will overlap, but more typically the fragments loaded in the machine are longer, and the bases in
the middle of the fragment are never read.

The sequence reads from the machine are saved into a file format named [FASTQ](
https://en.wikipedia.org/wiki/FASTQ_format) which contains both the sequence (of *ATCG* letters)
and the per-base quality score, which is an estimate of the error rate recorded by the machine as
it runs. Average quality drops off as the sequencer reads further into the fragment. It is
generally desirable to discard low-quality data, either by trimming off bad bases from the end or
discarding the whole bad read. There may also be "adapter sequences" within the reads which are
artefacts of preparing the DNA fragments for sequencing, not part of the biological sample. Tools
exist to detect and remove these from the FASTQ files.

In the most common case, the next step in DNA analysis after quality filtering is to map the
reads onto a known reference genome, which is the job of [an aligner](
https://en.wikipedia.org/wiki/List_of_sequence_alignment_software#Short-read_sequence_alignment).
This finds the most likely location where the raw read came from, allowing for a certain degree
of sequence mismatch. The second most common idea is to perform a [de-novo assembly](
https://en.wikipedia.org/wiki/De_novo_sequence_assemblers), putting the short reads together like
a jigsaw puzzle to try and make something as close to the complete chromosomes as possible. In
practise this process often produces very fragmentary results, yielding stretches of contiguous DNA
known as "contigs".

The sequencing of RNA (in practise this is [reverse transcribed](
https://en.wikipedia.org/wiki/Reverse_transcriptase) into DNA prior to being loaded into the
sequencer) gives an insight into gene activity. As with DNA, the analysis of these RNA reads also
involves aligning them to a reference genome, and in very simple terms we then count up how many
reads come from each gene in that genome. This shows which genes were active, and thus what
proteins were being produced, in the cells at the time the RNA was extracted. It's also possible
to perform de-novo assembly on RNA reads, attempting to reconstruct the sequences of the whole
genes and other RNA transcripts (aka. the transcriptome).

