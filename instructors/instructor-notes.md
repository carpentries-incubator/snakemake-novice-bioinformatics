---
title: Instructor Notes
---

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

## Questions from Learners

### Q: Example question here?

A: Example answer here.

