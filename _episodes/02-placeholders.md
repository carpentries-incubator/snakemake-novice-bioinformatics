---
title: "Placeholders and wildcards"
teaching: 0
exercises: 0
questions:
- "How do I make a generic rule?"
objectives:
- "Use Snakemake to count the lines in any file"
keypoints:
- "Add key points"
---

Similar to the workflows-snakemake thing that already exists. It intros wildcards after getting the first
workflow running.

In the previous ep you wrote two rules to count the sequences in two files. These work, but they are not
a very efficient use of Snakemake. We don't want to copy-paste rules avery time we want to process a new
input file.

{input} and {output} placeholders.
Use only in the 'shell' command.
Let's put them into the existing snakefile.

do it. test it.

Now add wildcards. You make up the names.
Use only in the 'input' and 'output' parts (actually there are other parts you can add but we'll come to them later)
A wildcard represents a variable part of the filename.
Let's collapse out rule into one rule.

do it. test it.

Now snakemake doesn't know what output we want, so we have to set it manually on the command line.

Maybe introduce the filtering here too? But not the chaining of rules yet.

Exercises:

Some multiple choice? Like what's wrong with this Snakefile?

What if I wanted to make the output be "gut_count_1.txt"? How could I make this work.



{% include links.md %}

