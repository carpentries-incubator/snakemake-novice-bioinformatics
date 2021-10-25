---
title: "Cleaning up"
teaching: 20
exercises: 10
questions:
- "How do I save disk space by removing temporary files?"
- "How do I protect important outputs from deletion?"
objectives:
- "Understand the function of *temporary* and *protected* outputs."
- "Learn about running in *touch* mode"
keypoints:
- "Cleaning up working files is good practise"
- "Make use of the `temporary()` function on outputs you don't need to keep"
- "Shadow rules can solve some issues with commands that produce unwanted files, but are not normally needed"
---
*For reference, [this is the Snakefile](../code/ep11.Snakefile) you should have to start the episode.
We're going back the original RNA-Seq workflow, not the assembly workflow from episode 11.*

## Temporary files in Snakemake

Analyses in bioinformatics inevitably generate a lot of files. Not all of these need to be kept, if they
are intermediate files produced during the workflow. Even if the files are not completely redundant, knowing
that you can regenerate them again with a single Snakemake command means you may well want to clean them up.
This way you save disk space and reduce clutter.

Remember that Snakemake only re-generates intermediate files if it needs them to make a target, so simply
deleting them is fine. Assuming you already made the *multiqc* report:

~~~
$ rm trimmed/ref*.fq
$ snakemake -j1 -p multiqc_out
...
Nothing to be done.
~~~

To get Snakemake to delete the files for you, mark them with the `temporary()` function. Much like the
`directory()` function, this is applied only on the outputs of rules. Any file marked as *temporary* will
be removed by Snakemake as soon as it is no longer needed.

To provide a better example, lets say we decide to compress the trimmed reads with *gzip*. It's very common
to store FASTQ files in this compressed format, and most software will read the files directly. Add a new
rule like so:

~~~
rule gzip_fastq:
    output: "{afile}.fq.gz"
    input:  "{afile}.fq"
    shell:
        "gzip -nc {input} > {output}"
~~~

Note that this will compress any `.fq` file in any subdirectory, because Snakemake wildcards can match full
paths. Now modify just the *salmon_quant* rule to work on the compressed files.

~~~
rule salmon_quant:
    output: directory("salmon.{sample}")
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
        fq1   = "trimmed/{sample}_1.fq.gz",
        fq2   = "trimmed/{sample}_2.fq.gz",
    ...
~~~

We'll pretend for now that Kallisto and FastQC will only work on the uncompressed files, even though these tools can
read `.fq.gz` files fine. Finally, declare the output of *trimreads* to be *temporary*.

~~~
rule trimreads:
  output: temporary("trimmed/{sample}.fq")
~~~

And now re-run the workflow. Since we modified the *trimreads* rule, we'll force that rule to be re-run with the
`-R` flag:

~~~
$ snakemake -j1 -p -Rtrimreads multiqc_out
...
$ ls trimmed/
etoh60_1_1.fq.gz  ref_1_1.fq.gz  temp33_1_1.fq.gz
etoh60_1_2.fq.gz  ref_1_2.fq.gz  temp33_1_2.fq.gz
etoh60_2_1.fq.gz  ref_2_1.fq.gz  temp33_2_1.fq.gz
etoh60_2_2.fq.gz  ref_2_2.fq.gz  temp33_2_2.fq.gz
etoh60_3_1.fq.gz  ref_3_1.fq.gz  temp33_3_1.fq.gz
etoh60_3_2.fq.gz  ref_3_2.fq.gz  temp33_3_2.fq.gz
~~~

Snakemake kept the uncompressed trimmed reads long enough to run *kallisto_quant* on all the samples, then removed
them leaving only the gzipped versions.

## Protected outputs and *touch* mode

One annoying thing about Snakemake is that, in moderately complex workflows, it may seem determined to re-run
a job for reasons that are unclear to you. For example, after adding the *gzip_fastq* rule above we re-ran
all of the *kallisto* and *salmon* quantifications, but with real data this could take hours or even days.

Let's pretend we just decided to *gzip* all the trimmed reads manually, ouside of Snakemake. This results
in files with new timestamps, so to simulate this we can just *touch* all the existing files.

~~~
$ touch trimmed/*
$ snakemake -n -p multiqc_out
...
Job counts:
	count	jobs
	1	multiqc
	9	salmon_quant
	10
~~~

Snakemake wants to re-run all the *salmon_quant* jobs, which makes sense, but we know the results are good, and
don't want to waste time re-making them, so we can fudge the timestamps using the `--touch` option. In the words
of the Snakemake manual, "This is used to pretend that the rules were executed, in order to fool future
invocations of snakemake."

~~~
$ snakemake -j 1 --touch -p multiqc_out
~~~

Another feature is the `protected()` function. This is rather like the opposite of the `temporary()` function
and says that once an output has been produced it must not be overwritten. In practise, Snakemake does this
by revoking write permissions on the files (as in `chmod -w {output}`).

> ## Exercise
>
> The longest operations in our test workflow are the genome indexing steps.
>
> Get Snakemake to protect the outputs from these steps. What happens if you try to overwite the protected files,
> by re-running the whole workflow with the `-F` option?
>
> > ## Answer
> >
> > ~~~
> > rule salmon_index:
> >     output:
> >         idx = protected(directory("{strain}.salmon_index"))
> >     ...
> >
> > rule kallisto_index:
> >     output:
> >         idx = protected("{strain}.kallisto_index"),
> >         log = "{strain}.kallisto_log",
> >     ...
> > ~~~
> >
> > When Snakemake is re-run, it will start processing the workflow and only fail when it comes to a protected file.
> > In the current version, Snakemake does not detect that the file is protected prior to starting executing jobs.
> > That may change in future versions.
> >
> {: .solution}
{: .challenge}

> ## Note
>
> In practise, it can be annoying to have Snakemake trying to regenerate files that have been protected so if you don't
> expect to be rebuilding indexes you could just remove or comment out the corresponding rules from the Snakefile,
> so Snakemake will not try to run them at all.
>
> Alternatively, once you have generated an important output, move the file away from your working directory. That
> way it will be easier not to accidentally clobber it (which is remarkably easy to do!).
>
{: .callout}

## Shadow rules

We'll not cover these in any detail but mention them for completeness. The
[Snakemake manual describes](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#shadow-rules)
the *shadow* directive very nicely.

In *shadow* mode, Snakemake links the input files into a temporary working directory, then runs the *shell*
command and finally copies the outputs back. If you have ever used [NextFlow](https://nextflow.io) this idea
will be familiar as NextFlow runs all workflow steps this way.

Advantages of using shadow rules are:

* Extra temporary files created by applications do not require explicit removal
* When running jobs in parallel (eg. `-j 2`) certain conflicts are resolved by not running multiple jobs in the same directory at once

Disadvantages are:

* Can be confusing when error messages reference the shadow directory
* Symlinks to subdirectories do not always work the way you expect
* Shadow directories are not always removed cleanly if Snakemake exits with an error

{% include links.md %}
