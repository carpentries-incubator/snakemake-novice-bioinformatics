---
title: "Globs and expansions"
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

Consider the following:

    for afile in reads/*.fq ; do
        outfile=$( basename $afile ).count
        wc -l $afile > $afile.count
    done



Note that we can also use Python syntax. But we don't expect Python knowledge.

{% include links.md %}

