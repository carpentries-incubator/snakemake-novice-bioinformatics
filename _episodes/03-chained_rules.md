---
title: "Chaining rules"
teaching: 0
exercises: 0
questions:
- "How do I combine rules into a workflow?"
- "How do I make a rule with multiple inputs and outputs?"
objectives:
- "Use Snakemake to filter and count the lines in a FASTQ file"
keypoints:
- "Add key points"
---

We now have a filter rule and a count rule. Consider how to combine them.

fastq_quality_trimmer -t 20 -l 100 -o filtered_etoh60_1_2.fq <reads/etoh60_1_2.fq

Chaining rules is simple if the names are chosen correctly. There's a bit of an art to it.
Demo is here.


Now we have the chaining concept. Yay. Can introduce the third alignment step here, and
a simple quality metric, so we have four rules by the end. Note that the alignment has multiple
inputs, so we have named inputs. Do I want to introduce this now or in placeholders part?
See where it fits.

Note that GG is running MultiQC and FastQC to we can also maybe bolt these on here?

{% include links.md %}

