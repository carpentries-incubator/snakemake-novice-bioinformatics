---
title: Cleaning up
teaching: 20
exercises: 15
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand the function of *temporary* outputs.
- Learn about running in *\--touch* mode
- Learn about *shadow* mode rules

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I save disk space by removing temporary files?
- How do I protect important outputs from deletion?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the final Snakefile from episodes 1 to 7](files/ep07.Snakefile) you may
use to start this episode.*

## Temporary files in Snakemake

Analyses in bioinformatics inevitably generate a lot of files. Not all of these need to be kept, if
they are intermediate files produced during the workflow. Even if the files are not completely
redundant, knowing that you can regenerate them again with a single Snakemake command means you may
well want to clean them up. This way you save disk space and reduce clutter.

Remember that Snakemake only re-generates intermediate files if it needs them to make a target,
so simply deleting intermediate files manually works fine. Assuming you already made the
*multiqc* report, and have not edited the rules in the Snakefile since:

```output
$ rm trimmed/ref*.fq
$ snakemake -j1 -p multiqc_out
...
Nothing to be done.
```

To get Snakemake to clean up files for you, mark them with the `temporary()` function. Much like
the `directory()` function, this is applied only on the outputs of rules. Any file marked as
*temporary* will be removed by Snakemake as soon as it is no longer needed.

To provide a better example, lets say we decide to compress the trimmed reads with *gzip*. It's
very common to store FASTQ files in this compressed format, and most software will read the gzipped
files directly, so there is no need to keep the uncompressed files. Add a new rule like so:

```source
rule gzip_fastq:
    output: "{afile}.fq.gz"
    input:  "{afile}.fq"
    shell:
        "gzip -nc {input} > {output}"
```

Note that this will compress any `.fq` file in any subdirectory, because Snakemake wildcards can
match full paths. Now modify just the *salmon\_quant* rule to work on the compressed files.

```source
rule salmon_quant:
    output: directory("salmon.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
        fq1   = "trimmed/{sample}_1.fq.gz",
        fq2   = "trimmed/{sample}_2.fq.gz",
    ...
```

Finally, declare the output of the *trimreads* rule to be *temporary*.

```source
rule trimreads:
    output: temporary("trimmed/{sample}.fq")
```

And now re-run the workflow. Since we modified the *trimreads* rule, Snakemake should see that it
needs to be rerun:

```output
$ snakemake -j1 -p multiqc_out
...
$ ls trimmed/
etoh60_1_1.fq.gz  ref_1_1.fq.gz  temp33_1_1.fq.gz
etoh60_1_2.fq.gz  ref_1_2.fq.gz  temp33_1_2.fq.gz
etoh60_2_1.fq.gz  ref_2_1.fq.gz  temp33_2_1.fq.gz
etoh60_2_2.fq.gz  ref_2_2.fq.gz  temp33_2_2.fq.gz
etoh60_3_1.fq.gz  ref_3_1.fq.gz  temp33_3_1.fq.gz
etoh60_3_2.fq.gz  ref_3_2.fq.gz  temp33_3_2.fq.gz
```

Snakemake kept the uncompressed trimmed reads long enough to run *kallisto\_quant* on all the
samples, then removed them leaving only the gzipped versions.

:::::::::::::::::::::::::::::::::::::::  challenge

## Working with gzipped files

Like Salmon, Kallisto can read compressed `.fq.gz` files directly. Amend the *kallisto\_quant*
rule to use gzipped input files.
Show that the new rule can be re-run on all samples without triggering the *trimreads* or
*gzip\_fastq* rules.

:::::::::::::::  solution

## Solution

In the first instance, this involves changing two input lines in the *kallisto\_quant*
rule definition.

```source
    fq1   = "trimmed/{sample}_1.fq.gz",
    fq2   = "trimmed/{sample}_2.fq.gz",
```

We can check the result by trying a *dry run* (`-n` option):

```bash
$ snakemake -n -p multiqc_out
```

You should see that Snakemake wants to run the *kallisto_quant* and *multiqc* steps but
no others.

::::::::::::::::::::::

## Removing the HTML reports from FastQC

We have no use for the HTML reports produced by FastQC. Modify the Snakefile to automatically
remove them.

:::::::::::::::  solution

## Solution

Amend the *html* output of the *fastqc* rule to mark it as `temporary()`:

```source
html = temporary("{indir}.{sample}_fastqc.html"),
```

Since the files are already there, Snakemake will not normally remove them unless the jobs
are re-run, so you could do that as was done for *trimreads* above.
However, there is also a *\--delete-temp-output* option which forces all temporary files in
the DAG to be removed, and this provides the cleanest way to remove the files after modifying
the Snakefile.

```bash
$ snakemake -p -j1 --delete-temp-output multiqc
```

::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Running in *\--touch* mode

One annoying thing about Snakemake is that, in moderately complex workflows, it may seem determined
to re-run a job for reasons that are unclear to you. For example, after adding the *gzip\_fastq*
rule above we re-ran all of the *kallisto* and *salmon* quantifications, but with real data this
could take hours or even days.

Let's pretend we just decided to *gzip* all the trimmed reads manually, ouside of Snakemake. This
results in files with new timestamps, so to simulate this we can just *touch* all the existing
files.

```bash
$ touch trimmed/*
$ snakemake -n -p multiqc_out
...
Job counts:
    count   jobs
    1       multiqc
    9       salmon_quant
    10
```

Snakemake wants to re-run all the *salmon\_quant* jobs (and the *kallisto\_quant* jobs, if you
completed the exercise above), which makes sense. However, we know the results are good, and don't
want to waste time re-making them, so we can fudge the timestamps using the *\--touch* option. In
the words of the Snakemake manual, "This is used to pretend that the rules were executed, in order
to fool future invocations of snakemake."

```bash
$ snakemake -j 1 --touch -p multiqc_out
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Protecting specific files

A related feature is the `protected()` function. This is rather like the opposite of the
`temporary()` function and says that once an output has been produced it must not be overwritten.
In practise, Snakemake does this by revoking write permissions on the files
(as in `chmod -w {output}`).

This works, but can be annoying because Snakemake will refuse to run if it believes it needs to
re-create a file which is protected. An alternative suggestion is, once you have generated
an important output file, move or copy the file away from your working directory. That way it will
be harder to accidentally clobber it (which is remarkably easy to do!).

::::::::::::::::::::::::::::::::::::::::::::::::::

## Shadow mode rules

Setting your rules to run in *shadow* mode helps clean up temporary files that are only used within
the job. The [Snakemake manual describes
](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules)
the `shadow` directive very nicely, so refer to that for more details.

Briefly, in *shadow* mode, Snakemake links the input files into a temporary working directory, then
runs the *shell* command and finally copies the outputs back. If you have ever used
[NextFlow](https://nextflow.io) this idea will be familiar as NextFlow runs all workflow steps
this way.

This test rule serves to demonstrate the operation of rules using the *shadow* directive.
The file `temp_file.txt` will not be there after the job has run, but the `shadow_out.txt` file
will be there because Snakemake sees that it is an output file and moves it back to the real
working directory,

```source
rule shallow_rule:
    output: "shadow_out.txt"
    shadow: "minimal"
    shell:
        """echo minimal shadow mode
           touch shadow_out.txt
           touch temp_file.txt
           tree `pwd`
        """
```

Advantages of using shadow rules are:

- Any temporary files created by applications do not require explicit removal
  (as with `temp_file.txt` in the above example).
- When running jobs in parallel (eg. `-j2`) certain conflicts related to temporary files will be
  avoided by not ever running multiple jobs in the same directory at once - we've already seen
  a case back in [Episode 06](06-awkward_programs.html) where the final version of the *fastqc*
  rule has this problem.

Disadvantages are:

- It can be confusing when error messages reference the shadow directory.
- Symlinks to subdirectories do not always work the way you expect.
- Shadow directories are not always removed cleanly if Snakemake exits with an error.

You may want to test your rules in normal mode first, then add `shadow: "minimal"` before you run
the workflow for real.


::::::::::::::::::::::::::::::::::::::::  challenge

## Removing the HTML reports (again)

Amend the *fastqc* rule once more so that the HTML reports are not even mentioned in the rule, and
will not appear in the working directory.

:::::::::::::::  solution

## Solution

The `shadow: "minimal"` directive will do the job nicely. You also need to remove mention of the
`.html` file from the list of outputs and the shell commands.

```code
rule fastqc:
    shadow: "minimal"
    output:
        zip  = "{indir}.{myfile}_fastqc.zip"
    input:  "{indir}/{myfile}.fq"
    shell:
        """fastqc -o . {input}
           mv {wildcards.myfile}_fastqc.zip {output.zip}
        """
```

In this case, marking the *html* output as `temporary()` or simply removing the file within the
*shell* part does work fine, but the good thing about the *shadow* approach is you do
not need to deal with or mention the unwanted file at all.

:::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is a Snakefile](files/ep12.Snakefile) incorporating the changes made in
this episode.*


:::::::::::::::::::::::::::::::::::::::: keypoints

- Cleaning up working files is good practise
- Make use of the `temporary()` function on outputs you don't need to keep
- Shadow rules can solve issues with commands that produce unwanted files

::::::::::::::::::::::::::::::::::::::::::::::::::

