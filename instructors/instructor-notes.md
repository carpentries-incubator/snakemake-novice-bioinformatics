---
title: Instructor Notes
---

# General

## How to introduce Snakemake

Prior to beginning the first lesson you want to say something about Snakemake in general.
As an instructor, saying how you yourself came across Snakemake and how you use it in your own
work is probably the best approach.

Otherwise, the info on [https://snakemake.readthedocs.io](https://snakemake.readthedocs.io) should
have everything you need, and the [rolling paper](https://f1000research.com/articles/10-33/v2)
has a nice graphic (fig. 2) showing the history of the Snakemake project. As of July 2024 this
paper has over 1000 citations.

You can also show the workflows available in the public catalogues to give an overview of what
people are doing with Snakemake.

* [Snakemake workflows in WokflowHub](
   https://workflowhub.eu/workflows?filter%5Bworkflow_type%5D=snakemake&order=downloads_desc)
* [Snakemake workflow catalogue](
   https://snakemake.github.io/snakemake-workflow-catalog/)

## When to use a workflow system?

Learners may ask when it is appropriate to use a system like Snakemake.

The paper [Workflow systems turn raw data into scientific knowledge](
https://pubmed.ncbi.nlm.nih.gov/31477884/) has a view on this:

::::::::::::::::::::::::::::::::::::::  discussion

So, do you need a workflow system? Not every task requires one, and there is a learning curve.
Scripting usually suffices for one-off tasks and when working out the pipeline itself. The
tipping point, most agree, comes when you need to run the same workflow over and over again, or
if the data are likely to be published.


::::::::::::::::::::::::::::::::::::::::::::::::::

I'd tend to disagree, in that some simple tasks work very nicely with Snakemake, and are easily
defined in one or two rules. Once you understand the fundamentals you are likely to start using
Snakemake for even these simple tasks.

Having said this, not every data analysis task is suited to Snakemake, or in some cases you may
only want to use Snakemake for part of a task, and do the rest with, say, regular scripting.

## Which is the best workflow system to use?

Snakemake! 🐍

But, in seriousness, other workflow systems are available. Some are better suited to different
tasks, and some users have a preference for one over another. For a large task, it is worth
investigating multiple options before committing to an approach.
[This GIT repository and associated paper](https://github.com/GoekeLab/bioinformatics-workflows)
comparing eight workflow systems is a good place to start. And in fact the previously mentioned
[Snakemake rolling paper](https://f1000research.com/articles/10-33/v2) compares Snakemake to
several other workflow systems.

You should also look through existing workflows on resources like [WorkflowHub](
https://workflowhub.eu), as someone may have already solved all or part of your problem.

## About the sample data files

The sample data files are not biologically meaningful. They are not even really from three
different conditions, so do not try to do any sort of real data analysis on them. For the purposes
of this course, the fact that they align to the genome is sufficient to make them a reasonable toy
dataset. At some point we could consider replacing this dataset with a new one, but it wouldn't
really add anything to the course.

It's possible that a learner will accidentally delete or overwite the input files. In this case,
note that a copy is available to download - see the link on [the setup page](../learners/setup.md).

## Choice of bioinformatics software

Like the toy dataset, the tools in this course are chosen to illustrate the workings of Snakemake.
The choice of older and simpler tools like *fastx toolkit* is deliberate, and reduces the burden of
maintenance of this course material as tools are updated.

In practise, learners may ask to go into more depth on the choice, configuration, and functionality
of the bioinformatics software. If you have the time and are confident talking about this then do
so, but if not then it is valid to reiterate that the focus of the course is on the orchestration
of analysis steps with Snakemake, not the choice of what software is best for any given analysis.

# Notes on specific episodes

## Episode 01 - Running commands with Snakemake

### Use of the -F flag

In the first few episodes we always run Snakemake with the `-F` flag, and it's not explained what
this does until Ep. 04. The rationale is that the default Snakemake behaviour when pruning the DAG
leads to learners seeing different output (typically the message "nothing to be done") when
repeating the exact same command. This can seem strange to learners who are used to scripting and
imperative programming.

The internal rules used by Snakemake to determine which jobs in the DAG are to be run, and which
skipped, are pretty complex, but the behaviour seen under `-F` is much more simple and consistent;
Snakemake simply runs every job in the DAG every time. You can think of `-F` as disabling the lazy
evaluation feature of Snakemake, until we are ready to properly introduce and understand it.

### Use of redirection (<) and shell arithmetic

A command is presented to count the sequences in a FASTQ file:

```
$ echo $(( $(wc -l <file.fq) / 4 ))
```

Understanding this in depth involves some advanced shell concepts that learners will not
necessarily be familiar with. However, other alternatives involve extra software, or use curly
brackets (which would have to be doubled-up) or are not robust (eg. `grep '^@' file.fq`).

## Episode 03 - Chaining rules

### Illustrating the wildcard matching process

There is a figure to illustrate the way Snakemake finds rules by wildcard matching and then tracks
back until it runs out of rule matches and finds a file that it already has. You may find that
showing an animated version of this is helpful, in which case
[there are some slides here](
https://github.com/carpentries-incubator/snakemake-novice-bioinformatics/files/9299078/wildcard_demo.pptx).

### Named inputs versus lists of inputs

In this course, we introduce named inputs and outputs before lists of inputs and outputs. This
results in shell commands like:

`"kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"`

Rather than the less readable version with a simple list of inputs:

`"kallisto quant -i {input[0]} -o kallisto.{wildcards.sample} {input[1]} {input[2]}"`

Later, we introduce lists of inputs in tandem with the `expand()` function. Of course it is
possible to have a list of outputs, but this is uncommon and not needed to solve any of the
challenges in this course. In fact, introducing lists of outputs may confuse learners as they
may think it is possible for a rule to yield a variable number of outputs in the manner of the old
`dynamic()` behaviour, which is not a thing.

## Episode 07 - Finishing the basic workflow

The episode describes some tactics for incorporating "awkward" programs and then mentions wrappers
as an aside at the end.

Wrappers are great when they work, and potentially infuriating when they do not. Also, users of
Snakemake are liable to come across tools that are not in the wrappers repository, or they may
even aim to contribute to this effort, in which case they need to understand the principle of
what is going on inside.

Converting the workflow to use wrappers was mostly straightforward, but here are some caveats:

1) The wrappers are designed for the specific versions of the tools specified in the Conda
requirements. It took some trial an error to find the right version of the Kallisto wrapper to
work with our older version of Kallisto (which was chosen for CPU compatibility).

2) All tool wrappers have sample code but it's not necessarily obvious what you may change
(normally wildcard names) and what you can't (input and output names).

3) Some wrappers use lists of inputs and outputs while others use named inputs and outputs. This
course has urged that outputs should always be named, but with wrappers the user must use
whatever setup the wrapper used.

## Episode 09 - Performance

### Running on cluster and cloud

Running workflows on HPC or Cloud systems could be a whole course in itself. The topic is too
important not to be mentioned here, but also complex to teach because you need a cluster to work
on.

The original author of this material would demonstrate running the workflow on the [Cirrus](
https://cirrus.ac.uk) system in Edinburgh. If you are teaching this lesson and have institutional
HPC then ideally you should liaise with the administrators of the system to make a suitable
installation of a recent Snakemake version and a profile to run jobs on the cluster job scheduler.
In practise this may be easier said than done!

If you are able to demonstrate Snakemake running on cloud as part of one of these courses then
we'd much appreciate any feedback on how you did this and how it went.


