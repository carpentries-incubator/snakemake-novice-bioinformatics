---
permalink: index.html
site: sandpaper::sandpaper_site
---

A lesson introducing the [Snakemake workflow system](https://snakemake.github.io/) for
bioinformatics analysis. Snakemake enables the writing of reliable, scalable and
reproducible scientific workflows as a series of chained rules. Simple workflows to replace shell
scripts can be written in a few lines of code, while for more complex cases there is support for
conda integration, software containers, cluster execution, cloud execution, etc. You can also add
Python and R code directly into your workflow.

This lesson introduces the core concepts of Snakemake in the context of a typical analysis task,
aligning short cDNA reads to a reference transcriptome. Later episodes focus on practical questions
of workflow design, debugging and configuration.

We also look at the Conda integration feature of Snakemake, with which you can author reproducible
and shareable workflows with a fully-specified software environment.

In the planning phase of writing this course material we outlined some [learner profiles
](profiles.html), to expand on who we think will benefit from this lesson and why.

<!-- this is an html comment -->

[comment]: # (This is a markdown comment and will not be rendered into the HTML at all)

::::::::::::::::::::::::::::::::::::::::::  prereq

## Learner Prerequisites

See the [prerequisites](prereqs.html) page for a full list of skills and concepts we assume that
learners will know prior to taking this lesson. In brief:

This is an intermediate lesson and assumes learners have some prior experience in bioinformatics:

- Familiarity with the [Bash command shell](https://swcarpentry.github.io/shell-novice), including
  concepts like pipes, variables, loops and scripts.
- Knowing about bioinformatics fundamentals like the [FASTQ file format
  ](https://en.wikipedia.org/wiki/FASTQ_format) and [read mapping
  ](https://en.wikipedia.org/wiki/Read_\(biology\)#NGS_and_read_mapping),
  in order to understand the example workflows.

No previous knowledge of Snakemake or workflow systems, or Python programming, is assumed.


::::::::::::::::::::::::::::::::::::::::::::::::::


