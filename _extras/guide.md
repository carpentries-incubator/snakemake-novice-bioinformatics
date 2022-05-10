---
title: "Instructor Notes"
---

# General

## How to introduce Snakemake

Prior to beginning the first lesson you want to say something about Snakemake. As an instructor,
saying how you yourself came across Snakemake and how you use it in your own work is probably
the best approach.

Otherwise, the info on https://snakemake.readthedocs.io should have everything you need, and the
"rolling paper" has a nice graphic showing the history of the Snakemake project.

## When to use a workflow system?

Learners may ask when it is appropriate to use a system like Snakemake.

The paper [Workflow systems turn raw data into scientific knowledge](https://pubmed.ncbi.nlm.nih.gov/31477884/) has a view on this:

> So, do you need a workflow system? Not every task requires one, and there is a learning curve. Scripting
> usually suffices for one-off tasks and when working out the pipeline itself. The tipping point, most agree,
> comes when you need to run the same workflow over and over again, or if the data are likely to be published.

I'd tend to disagree, in that some simple tasks work very nicely with Snakemake, and are easily defined in one
or two rules. Once you understand the fundamentals you are likely to start using Snakemake for even these simple
tasks.

Having said this, not every data analysis task is suited to Snakemake, or in some cases you may only want to use
Snakemake for part of a task, and do the rest with regular scripting.

## Which is the best workflow system to use?

Snakemake.

But, in seriousness, other workflow systems are available. Some are better suited to different tasks, and some users
have a preference for one over another. For a large task, it is worth investigating multiple options before
committing to an approach.

# Notes on specific episodes

## TODO

{% include links.md %}
