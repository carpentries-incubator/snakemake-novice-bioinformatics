---
title: How Snakemake plans what jobs to run
teaching: 40
exercises: 30
---

::::::::::::::::::::::::::::::::::::::: objectives

- View the DAG for our pipeline
- Understand the logic Snakemake uses when running and re-running jobs

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I visualise a Snakemake workflow?
- How does Snakemake avoid unecessary work?
- How do I control what steps will be run?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the Snakefile](files/ep03.Snakefile) you should have to start the
episode.*

## The DAG

You may have noticed that one of the messages Snakemake always prints is:

```output
Building DAG of jobs...
```

A DAG is a **Directed Acyclic Graph** and it can be pictured like so:

![][fig-dag]

The above DAG is based on our four existing rules, and shows all the jobs Snakemake would run to
trim, count and quantify the *ref1* sample.

:::::::::::::::::::::::::::::::::::::::  checklist

## Note that:

- A rule can appear more than once, with different wildcards (a **rule** plus **wildcard values**
  defines a **job**)
- The arrows show dependency ordering between jobs
- Snakemake can run the jobs in any order that doesn't break dependency - for example
  `kallisto quant` cannot run until
  both `kallisto index` and `trimreads` have completed, but it may run before or after
  `countreads`
- This is a work list, *not a flowchart*, so there are no if/else decisions or loops - Snakemake
  runs every job in the DAG exactly once
- The DAG depends both on the Snakefile *and* on the requested target outputs, and the files
  already present
- When building the DAG, Snakemake does not look at the *shell* part of the rules at all - only
  when running the DAG will Snakemake check that the shell commands are working and producing the
  expected output files

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## How many jobs?

If we asked Snakemake to run *kallisto\_quant* on all three of the reference samples
(ref1, ref2, ref3), how many jobs would that be in total?

:::::::::::::::  solution

## Solution

10 in total: 3 \* kallisto\_quant // 6 \* trim // 1 \* kallisto\_index // 0 \* countreads

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Snakemake is lazy, and laziness is good

For the last few episodes, we've told you to run Snakemake like this:

```bash
$ snakemake -j1 -F -p desired_output_file
```

As a reminder, the `-j1` flag tells Snakemake to run one job at a time, and `-p` is to print out
the shell commands before running them.

The `-F` flag turns on `forceall` mode, and in normal usage you don't want this.

At the end of the last chapter, we generated some kallisto results by running:

```bash
$ snakemake -j1 -F -p kallisto.temp33_1/abundance.h5
```

Now try without the `-F` option. Assuming that the output files are already created, you'll see
this:

```bash
$ snakemake -j1 -p kallisto.temp33_1/abundance.h5
Building DAG of jobs...
Nothing to be done.
Complete log: /home/zenmaster/data/yeast/.snakemake/log/2021-04-23T172441.519500.snakemake.log
```

In normal operation, Snakemake only runs a job if:

1. A target file you explicitly requested to make is missing
2. An intermediate file is missing and it is needed in the process of making a target file
3. Snakemake can see an input file which is newer than an output file

Let's demonstrate each of these in turn, by altering some files and re-running Snakemake without
the `-F` option.

```bash
$ rm kallisto.temp33_1/*
$ snakemake -j1 -p kallisto.temp33_1/abundance.h5
```

This just re-runs the kallisto quantification - the final step.

```bash
$ rm trimmed/temp33_*.fq
$ snakemake -j1 -p kallisto.temp33_1/abundance.h5
```

"Nothing to be done" - some intermediate output is missing but Snakemake already has the file you
are telling it to make, so it doesn't worry.

```bash
$ touch transcriptome/*.fa.gz
$ snakemake -j1 -p kallisto.temp33_1/abundance.h5
```

The `touch` command is a standard Linux command which sets the timestamp of the file, so now the
transcriptome looks to Snakemake as if it was just modified.

Snakemake sees that one of the input files used in the process of producing
`kallisto.temp33_1/abundance.h5` is newer than the existing output file, so it needs to run the
`kallisto index` and `kallisto quant` steps again. Of course, the `kallisto quant` step needs the
trimmed reads which we deleted earlier, so now the trimming step is re-run also.

## Explicitly telling Snakemake what to re-run

The default timestamp-based logic is really useful when you want to:

1. Change or add some inputs to an existing analysis without re-processing everything
2. Continue running a workflow that failed part-way

But it doesn't help us in the situation when rules in the Snakefile change, rather than input
files. Snakemake won't see that the results are out-of-date. For example, if we changed the quality
cutoffs within the trimreads rule then Snakemake would not automatically re-run those rules,
because it only checks that the output file is newer than the input file.

The `-R` flag allows you to explicitly tell Snakemake that a rule has changed and that all outputs
from that rule need to be re-evaluated.

```bash
$ snakemake -j1 -R trimreads -p kallisto.temp33_1/abundance.h5
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Note on `-R`

Due to a quirk of the way Snakemake parses command-line options, you need to make sure there are
options after the `-R ...`, before the list of target outputs. If you don't do this, Snakemake
will think that the target files are instead items to add to the `-R` list, and then when
building the DAG it will just try to run the default rule.

The easiest way is to put the `-p` flag before the target outputs. Then you can list multiple
rules to re-run, and also multiple targets, and Snakemake can tell which is which.

```bash
$ snakemake -j1 -R trimreads kallisto_index -p kallisto.temp33_1/abundance.h5 kallisto.temp33_2/abundance.h5
```

The reason for using the `-p` flag here is because you generally always want this option.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Another note on `-R`

In [Snakemake version 7.8
](https://github.com/snakemake/snakemake/blob/main/CHANGELOG.md#780-2022-05-24)
the behaviour is altered to automatically detect changes to rules, as well as parameters and
input file contents. This to some extent makes the `-R` flag redundant, but the exact behaviour
may change again in future as there are [problems with the logic in more complex cases
](https://github.com/snakemake/snakemake/issues/1694).
You should check the manual for whatever version of Snakemake you are using, and particularly the
`--rerun-triggers` option, to see what behaviour is expected.

Snakemake is under active development, so changes to behaviour and features are always liable to
be introduced in new versions.

::::::::::::::::::::::::::::::::::::::::::::::::::

The `-f` flag specifies that the target outputs named on the command line should always be
regenerated, so you can use this to explicitly re-make specific files.

```bash
$ snakemake -j1 -f -p kallisto.temp33_1/abundance.h5
```

This always re-runs *kallisto\_quant*, regardless of whether the output file is there already. For
all intermediate outputs, Snakemake applies the default timestamp-based logic. Contrast with `-F`
which runs the entire DAG every time.

## Visualising the DAG

Snakemake can draw a picture of the DAG for you, if you run it like this:

```bash
$ snakemake -f --dag kallisto.etoh60_1/abundance.h5 | gm display -
```

Using the `--dag` option implicitly activates the `-n` (dry-run) option so that Snakemake will
not actually run any jobs, it will just print the DAG and stop. Snakemake prints the DAG in a text
format so we use the `gm` command to make this into a picture and show it on the screen.

:::::::::::::::::::::::::::::::::::::::::  callout

## Note on `gm display`

The `gm` command is provided by the [GraphicsMagick toolkit](https://www.graphicsmagick.org/).
On systems where `gm` will not display an image directly, you may instead save it to a PNG file.
You will need the `dot` program from the [GraphViz package](https://graphviz.org/) installed.

```bash
$ snakemake -f --dag kallisto.etoh60_1/abundance.h5 | dot -Tpng > dag.png
```

::::::::::::::::::::::::::::::::::::::::::::::::::

![][fig-dag2]

The boxes drawn with dotted lines indicate steps that are not to be run, as the output files are
already present and newer than the input files.

:::::::::::::::::::::::::::::::::::::::  challenge

## Updating a parameter and re-running

Run `kallisto_quant` on all three of the **etoh60** samples. Now edit the Snakefile so that the
`qual_threshold` for `trimreads` is "22", rather than "20".

How would you get Snakemake to update all three Kallisto results:

1) By using the -R flag

2) By using the -f flag

3) By using the touch command

4) By deleting some of the existing files

Use the `--dag` option as shown above to check your answers.

:::::::::::::::::::: solution

## Solution

This will make all the Kallisto results in the first place:

```bash
$ snakemake -j1 -p kallisto.etoh60_{1,2,3}/abundance.h5
```

*The {1,2,3} syntax is expanded by the shell into the 3 file names. You could also type all
three names in full.*

1) -R flag
```bash
$ snakemake -Rtrimreads --dag kallisto.etoh60_{1,2,3}/abundance.h5 | gm display -
```

2) -f flag
```bash
$ snakemake -j1 -p --dag -f trimmed/etoh60_{1,2,3}_{1,2}.fq kallisto.etoh60_{1,2,3}/abundance.h5 | gm display -
```

3) touch command
```bash
$ touch reads/etoh60_*.fq`
```

4) by deleting files
```bash
$ rm -r trimmed/etoh60_*.fq kallisto.etoh60_*
```

:::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Removing files to trigger reprocessing

In general, getting Snakemake to re-run things by removing files is a bad idea, because it's easy
to forget about intermediate files that actually contain stale results and need to be updated.
Using the `-R` flag or `touch` command is simpler and more reliable. If in doubt, and if it will
not be too time consuming, keep it simple and just use `-F` to run the whole workflow from scratch.

For the opposite case where you want to avoid re-running particular steps, see the `‑‑touch`
option of Snakemake mentioned [later in the course.
](12-cleaning_up.md)

::::::::::::::::::::::::::::::::::::::::::::::::::



[fig-dag]: fig/dag_1.svg {alt='Diagram showing jobs as coloured boxes joined by arrows representing
data flow. The box labelled as kallisto\_index is in green at the top, with two blue boxes labelled
trimreads and two yellow boxes labelled countreads. The blue trimreads boxes have arrows into the
respective yellow countreads boxes. Finally there is a kallisto\_quant job shown as a red box, with
incoming arrows from both the trimreads box as well as the kallisto\_index box.'}
[fig-dag2]: fig/dag_2.png {alt='A DAG for the partial workflow with four boxes, representing two
trimreads jobs and a kallisto\_index job, then a kallisto\_quant job receiving input from the
previous three, The boxes for the kallisto\_index and trimreads jobs are dotted, but the
kallisto\_quant box is solid.'}


:::::::::::::::::::::::::::::::::::::::: keypoints

- A **job** in Snakemake is a **rule** plus **wildcard values** (determined by working back from the requested output)
- Snakemake plans its work by arranging all the jobs into a **DAG** (directed acyclic graph)
- If output files already exist, Snakemake can skip parts of the DAG
- Snakemake compares file timestamps to determine if outputs need regenerating

::::::::::::::::::::::::::::::::::::::::::::::::::


