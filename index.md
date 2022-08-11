---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---

A lesson introducing the [Snakemake workflow system](https://snakemake.github.io/) for bioinformatics
analysis. The Snakemake system enables the writing of reliable, scalable and reproducible scientific
workflows as a series of chained rules. Simple workflows to replace shell scripts can be written in a
few lines of code, while for more complex cases there is support for conda integration, software containers,
cluster execution, cloud execution, etc. You can also add Python and R code directly into your workflow.

This lesson introduces the core concepts of Snakemake in the context of a typical analysis task, aligning
short cDNA reads to a reference transcriptome. Later episodes focus on practical questions of workflow design,
debugging and configuration.

We also look at the Conda integration feature of Snakemake, with which you can author reproducible
and shareable workflows with a fully-specified software environment.

In the planning phase of writing this course material we outlined some [learner profiles](learner_profiles/),
to expand on who we think will benefit from this lesson and why.

<!-- this is an html comment -->

{% comment %} This is a comment in Liquid {% endcomment %}

> ## Prerequisites
>
> This is an intermediate lesson and assumes learners have some prior experience in bioinformatics:
> * Familiarity with the [Bash command shell](http://swcarpentry.github.io/shell-novice), including concepts
    like pipes, variables and loops.
> * Knowledge of bioinformatics fundamentals like the FASTQ file format and short read mapping,
>   in order to understand the example workflow.
>
> No previous knowledge of Snakemake or workflow systems, or Python programming, is required.
{: .prereq}

{% include links.md %}
