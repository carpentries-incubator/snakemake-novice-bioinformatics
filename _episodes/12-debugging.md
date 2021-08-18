---
title: "Planning and debugging your own workflows"
teaching: 0
exercises: 0
questions:
- "How do I approach making my own workflow?"
- "How do I visualize what Snakemake is doing?"
- "How do I understand and debug the errors I see?"
objectives:
- "Know the three stages of execution in Snakemake"
- "Debug some common errors in each phase"
keypoints:
- "Add key points"
---

## When designing a Snakemake workflow

* Start out with pen and paper. Sketch how you expect the DAG to look.
* Work out the inputs and the final outputs for the stage.
* List any parameters or settings that might need to be adjusted.
* Work out which rules (if any) aggregate or split inputs (these may need input functions).
* Get the `input` and `output` parts of your Snakemake rules right before worrying about the `shell` sections.
  Remember that Snakemake resolves all the steps before actually running the `shell` commands so you should
  --dryrun the workflow before running it for real, or even before you put the `shell` commands in.
* Snakemake starts matching rules from the target file back to the input files, so you might prefer to work
  on your rules in this same order, starting with the last step in each stage.
* Choose the names for your *input*s, *output*s and *{wildcards}* carefully to make your Snakefile as readable
  as possible.

TODO - At this point we could introdue a new workflow like the assembly workflow from before? Do we have time to talk
about velvet assembler?

## Debugging a Snakemake workflow

As with all programming jobs, you will inevitably see errors the first time you try to run a new Snakefile,
because it's very easy to make typos. Review the hints from yesterday and remember that you need to establish
which phase of execution failed, then look at the most common error causes for that phase.  If you get syntax
errors, look out for unexpected or missing text colouring in GEdit which might hint at the problem.
You should also check your punctuation: brackets, quotes, colons, line indentation and commas all need to be right.
If you can't see the problem, ask a demonstrator for help.

TODO - I need to make some broken example Snakefiles and then get participants to fix them. This is tricky.

Easier to have an open-ended exercise and let the problems come up naturally?

