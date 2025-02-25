---
title: Complex outputs, logs and errors
teaching: 30
exercises: 20
---

::::::::::::::::::::::::::::::::::::::: objectives

- Add an RNA quantification step in the data analysis
- Learn about adding log outputs to rules
- Understand why and how Snakemake deals with missing outputs

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How can we start to analyse the sample data?
- What can cause a job in Snakemake to fail?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the Snakefile](files/ep03.Snakefile) you should have to start
the episode.*

## Adding a transcript counting step to the pipeline

::::::::::::::::::::::::::::: instructor

The main point of the following exercise is to get some insight into what the learners might want
to do with Snakemake. You should focus on the answers received rather than dwelling on the sample
answers below.

::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Thinking about your own workflows

Think about any data processing task you have done yourself, and write down three or four steps
from that workflow.

What were the inputs to, and outputs from, each step?

How did the steps connect up, in terms of data going from one to the next? You may want to sketch
this out and use arrows to indicate the linkages between the steps.

:::::::::::::::  solution

## Solution

### A simple bioinformatics workflow

```
"Align reads to genome" : input=[reads, reference], output=BAM file

"Sort BAM file" : input=BAM file, output=sorted BAM file

"Summarize coverage" : input=[sorted BAM file, reference], output=histogram
```

### A workflow for a robot to brew a mug of tea

```
"Boil Water" : input=cold water, output=hot water

"Brew Tea" : input=[hot water, tea bag, mug], output=tea infusion

"Add Milk And Sugar" : input=[tea infusion, milk, sugar], output=cuppa
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Introducing Kallisto

Let's add another rule to our Snakefile. The reads we have are from a yeast RNA-seq experiment so
we might reasonably want to quantify transcript abundance using the **kallisto** program. The
command to do so looks like this:

```bash
$ kallisto quant -i index_file -o output_dir in_1.fastq in_2.fastq
```

This command has three input files:

1. The transcriptome index
2. The first of the paired FASTQ files
3. The second of the paired FASTQ files

And it produces a directory of output files. According to
[the Kallisto manual](https://pachterlab.github.io/kallisto/manual#quant) this directory will
have three output files in it:

1. abundance.h5
2. abundance.tsv
3. run\_info.json

We'll not worry about what the contents of these files mean just now, or how Kallisto generates
them. We just know that we want to run the `kallisto quant` command and have it make the output
files, and the output files are going to be useful once we add later steps in the analysis.

Making a rule with multiple inputs and outputs works much like we previously saw.

```source
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"
```

::::::::::::::::::::::::::::: instructor

## Named inputs versus lists of inputs

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

::::::::::::::::::::::::::::::::::

There are many things to note here:

 1. The individual input and output files are all given names.
 2. We've used the wildcard name `{sample}` rather than `{myfile}` because this will match only the
    sample name, eg `ref1`, not `ref1_1`. Snakemake doesn't care what name we use, but carefully
    chosen names make for more readable rules.
 3. Because `kallisto quant` only takes the output directory name, we've used the placeholder
    `{wildcards.sample}` rather than `{output}` which would expand to the full file names.
 4. We've chosen to only quantify the *trimmed* version of the reads.
 5. We don't actually have the `{input.index}` file yet. This will need to be created using the
    `kallisto index` command.

Even though the rule is not going to work without the index, we can still run it to check that
Snakemake is happy with the rule definition.

:::::::::::::::::::::::::::::::::::::::::  callout

## Running Kallisto on all replicates

If you are previously familiar with the Kallisto software, you may be thinking about running
Kallisto on all replicates of the condition at once. We'll look at this later in the course, but
for now we will be running Kallisto once per sample, ie. once for each pair of FASTQ files.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Running the kallisto\_quant rule

Given that the *index* input is missing, what would you expect Snakemake to do if the new rule
was run now?

Try it by telling Snakemake to run the new rule on the files `ref1_1.fq` and `ref1_2.fq`.
Since the rule defines multiple outputs, asking for any one of the output files will be enough.

::::::::::::::::::::::::::  solution

## Solution

```bash
$ snakemake -j1 -F -p kallisto.ref1/abundance.h5
```

Resulting error is "Missing input files". Note that Snakemake does not try to run kallisto at
all. It stops while still making a work plan, because it sees that a required input file is
missing.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Building the index

Instruct Snakemake how to build the genome index as part of the pipeline by adding another rule.
The command we need to run is:

```bash
$ kallisto index -i index_file_to_make fasta_file_to_index
```

The file to be indexed is `transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz`. As
there is only one input to the rule you don't have to give it a name, but you may do so if you
prefer.

Make it so that the terminal messages printed by the program are captured to a file, and therefore
your rule will have two separate outputs: the *index file* and the *messages file*. Note that the
program prints messages on *stderr*, so you will need to use `>&` rather than `>` to capture that
output.

:::::::::::::::  solution

## Solution

The rule you write could look something like this, but there are many variations that will work
just as well. Since there is only one transcriptome in the project, you may feel that use of
the `{strain}` wildcard is overkill, but who's to say we might not want to use another
in future?

```source
rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
        messages = "{strain}.kallisto_stderr",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {output.messages}"
```

::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::::

## Log outputs in Snakemake

All being well, our new rules are now ready to run Kallisto, and we can analyse any sample we
like.

```bash
$ snakemake -j1 -F -p kallisto.ref3/abundance.h5
...lots of output...
4 of 4 steps (100%) done
Complete log: /home/zenmaster/data/yeast/.snakemake/log/2021-04-23T142649.632834.snakemake.log
```

There is one particular improvement we can make, since Snakemake has a dedicated rule field for
outputs that are [log files](
https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files).
These are mostly treated the same as regular outputs except that log files are always kept even
if the job produces an error, so you can look at the log to help diagnose the error.

For an output to be treated as a log file, list it under `log:` instead of `output:` and then
within the shell command use the placeholder `{log}` instead of `{output}`.

::::::::::::::::::::::::::::::::::  challenge

## Using an explicit log output

Modify the solution to the previous challenge so that it uses the `log` keyword to capture the
terminal output from `kallisto quant`.

:::::::::::::::  solution

## Solution

```source
rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    log:
        messages = "{strain}.kallisto_stderr",
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {log.messages}"
```

The order of the `log:`, `output:` and `input:` parts can be however you like. In this case since
there is now a single input file, a single output file, and a single log file, you may feel that
there is no point naming them all.

```source
rule kallisto_index:
    output: "{strain}.kallisto_index"
    input:  "transcriptome/{strain}.cdna.all.fa.gz"
    log:    "{strain}.kallisto_stderr"
    shell:
        "kallisto index -i {output} {input} >& {log}"
```

:::::::::::::::::::::

::::::::::::::::::::::::::::::::


## Dealing with a *missing files* error

We'll end the chapter by looking at a common problem that can arise if you mistype a file
name in a rule. Remember that we wrote the rule based on the expected output filenames given in the
Kallisto manual. In an older version of this manual there was a typo where the file names were
incorrectly given as `abundances.h5` and `abundances.tsv`, with the extra `s` on each.

It may seem silly to break the workflow when we just got it working, but it will be instructive,
so edit the Snakefile and change these names to the incorrect versions.

```source
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundances.h5",
        tsv  = "kallisto.{sample}/abundances.tsv",
        json = "kallisto.{sample}/run_info.json",
...
```

To keep things tidy, this time we'll manually remove the output directory.

```bash
$ rm -rvf kallisto.ref1
```

And re-run. Note that the output file name you'll need to use on the command line must match the
edited Snakefile, or you will get a `MissingRuleException`.

```output
$ snakemake -j1 -F -p kallisto.ref1/abundances.h5

...
kallisto quant -i Saccharomyces_cerevisiae.R64-1-1.kallisto_index -o kallisto.ref1 trimmed/ref1_2.fq trimmed/ref1_2.fq

[quant] fragment length distribution will be estimated from the data
...more Kallisto output...
[   em] the Expectation-Maximization algorithm ran for 265 rounds

Waiting at most 5 seconds for missing files.
MissingOutputException in line 24 of /home/zenmaster/data/yeast/Snakefile:
Job Missing files after 5 seconds:
kallisto.ref1/abundances.h5
kallisto.ref1/abundances.tsv
This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.
Job id: 0 completed successfully, but some output files are missing. 0
  File "/opt/python3.7/site-packages/snakemake/executors/__init__.py", line 583, in handle_job_success
  File "/opt/python3.7/site-packages/snakemake/executors/__init__.py", line 259, in handle_job_success
Removing output files of failed job kallisto_quant since they might be corrupted:
kallisto.ref1/run_info.json
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
Complete log: /home/zenmaster/data/yeast/.snakemake/log/2021-04-23T142649.632834.snakemake.log
```

There's a lot to take in here. Some of the messages are very informative. Some less so.

1. Snakemake did actually run kallisto, as evidenced by the output from kallisto that we see on the
  screen.
2. There is no obvious error message being reported by kallisto.
3. Snakemake complains some expected output files are missing: `kallisto.ref1/abundances.h5` and
  `kallisto.ref1/abundances.tsv`.
4. The third expected output file `kallisto.ref1/run_info.json` was found but has now been
  removed by Snakemake.
5. Snakemake suggest this might be due to "filesystem latency".

This last point is a red herring. "Filesystem latency" is not an issue here, and never will be
since we are not using a network filesystem. We know what the problem is, as we deliberately caused
it, but to diagnose an unexpected error like this we would investigate further by looking at the
`kallisto.ref1` subdirectory.

```bash
$ ls kallisto.ref1/
abundance.h5  abundance.tsv
```

Remember that Snakemake itself does not create any output files. It just runs the commands you put
in the `shell` sections, then checks to see if all the expected output files have appeared.

So if the file names created by kallisto are not exactly the same as in the Snakefile you will get
this error, and you will, in this case, find that some output files are present but others
(`run_info.json`, which was named correctly) have been cleaned up by Snakemake.

:::::::::::::::::::::::::::::::::::::::::  callout

## Errors are normal

Don't be disheartened if you see errors like the one above when first testing your new Snakemake
pipelines. There is a lot that can go wrong when writing a new workflow, and you'll normally need
several iterations to get things just right. One advantage of the Snakemake approach compared to
regular scripts is that Snakemake fails fast when there is a problem, rather than ploughing on
and potentially running junk calculations on partial or corrupted data. Another advantage is that
when a step fails we can safely resume from where we left off, as we'll see in the next episode.

::::::::::::::::::::::::::::::::::::::::::::::::::

Finally, edit the names in the Snakefile back to the correct version and re-run to confirm that all
is well. Assure yourself that that the rules are still generic by processing the *temp33\_1*
sample too:

```bash
$ snakemake -j1 -F -p kallisto.ref1/abundance.h5  kallisto.temp33_1/abundance.h5
```

*For reference, [this is a Snakefile](files/ep04.Snakefile) incorporating the changes made in
this episode.*


:::::::::::::::::::::::::::::::::::::::: keypoints

- Try out commands on test files before adding them to the workflow
- You can build up the workflow in the order that makes sense to you, but always test as you go
- Use log outputs to capture the messages printed by programs as they run
- If a shell command exits with an error, or does not yield an expected output then Snakemake will
  regard that as a failure and stop the workflow

::::::::::::::::::::::::::::::::::::::::::::::::::


