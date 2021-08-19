---
title: "Handling awkward programs"
teaching: 20
exercises: 30
questions:
- "How do I handle real bioinformatics tools?"
- "How do I define a rule where the output is a directory?"
objectives:
- "Understand the different choices available when defining a rule"
- "Learn about the `directory()` function"
- "Add multiple commands to the shell section of a rule"
keypoints:
- "Add key points"
---

We've now seen how to link rules in a pipeline and how to merge results at the final step. This is the basic
pattern for many analysis workflows. For simplicity, we used mostly basic shell commands, but we'll now replace
these with bioinformatics tools:

* **Salmon** is a alternative to Kallisto, using a different alignment algorithm
* **FastQC** calculates a variety of metrics on a FASTQ file and produces an HTML report and a ZIP file.
* **MultiQC** combines the reports from various tools, including FastQC, Kallisto, and Salmon, into a single HTML report.

## The full workflow we are constructing

![Our full QC workflow][fig-workflow]

Real programs like this can have quirks like:

* Not allowing you to specify input or output file locations
* Outputting a whole directory of files
* Creating temporary files
* Not supporting robust error handling

With a little care we can handle all these problems.

## Adding a FastQC rule

We'll first run FastQC on a single input file and see what is produced.

~~~
$ fastqc reads/ref_1_1.fq
$ ls -ltr reads
...
-rw-r--r-- 1 zenmaster  users   464852 Jun  9 14:31 ref_1_1_fastqc.zip
-rw-r--r-- 1 zenmaster  users   654810 Jun  9 14:31 ref_1_1_fastqc.html
~~~

It's possible to supply multiple input files at once, but the resulting output is exactly the same as if processing the
files one at a time. Two files are produced for each FASTQ file, and these files appear in the same directory as the
input file.
The `fastqc` command does not let us specify the output filenames, but we can set the output directory name.

~~~
$ fastqc --help
...
   -o --outdir     Create all output files in the specified output directory.
                   Please note that this directory must exist as the program
                   will not create it.  If this option is not set then the
                   output file for each sequence file is created in the same
                   directory as the sequence file which was processed.
...
~~~

For the `countreads` rule we wrote earlier, we chose our preferred output file name first, then wrote the shell command
so as to put the results into that file. This allowed us to have a rule the counts both *trimmed* and *untrimmed* reads.

~~~
# Our existing countreads rule...
rule countreads:
  output: "{indir}.{asample}.fq.count"
  input:  "{indir}/{asample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~

To do the same with FastQC have four options:

1. Work with the default file names produced by FastQC and leave the reports in the same directory with the FASTQ files.
1. Make the outputs in a new directory named, eg. "reads.fastqc.ref_1_1/" (similar to what we did with Kallisto).
1. Do this, but explicitly tell Snakemake that the directory _is_ the output.
1. Force our preferred naming convention by renaming the FastQC output files within the rule.

We'll try all four, as all are valid options. It's often the case in Snakemake that you have options like this.

> ## Exercise - adding a FastQC rule using the default output file names
>
> Fill in the ??? to make a working rule for FastQC where `indir` may be "reads" or "trimmed". Do not change the
> shell command or input pattern at all. Remember FastQC always makes two output files, so add two named outputs.
>
> ~~~
> rule fastqc:
>   output:
>       ???
>   input:  "{indir}/{asample}.fq"
>   shell:
>       "fastqc {input}"
> ~~~
>
> > ## Solution
> >
> > Since the `shell` command is not to be changed, the output names will be dictated by fastqc.
> >
> > ~~~
> > rule fastqc:
> >   output:
> >       html = "{indir}/{asample}_fastqc.html",
> >       zip  = "{indir}/{asample}_fastqc.zip",
> >   input:  "{indir}/{asample}.fq"
> >   shell:
> >       "fastqc {input}"
> > ~~~
> {: .solution}
{: .challenge}

This rule is fine, but maybe we don't want to put the reports in with the sequences. As a general principle, when
writing Snakemake rules, you want to be in charge of output file names. FastQC lets you specify the output
directory, so we can use that.

> ## Exercise - a FastQC rule where the output files go into a new directory
>
> Modify the rule so that the output files go into a new directory. This will be very similar to the rule for
> `kallisto quant`.
>
> For example, when running on the file "trimmed/ref_1_1.fq" the outputs will be
> ~~~
> trimmed.fastqc.ref_1_1/ref_1_1_fastqc.html
> trimmed.fastqc.ref_1_1/ref_1_1_fastqc.zip
> ~~~
>
> > ## Solution
> >
> > This involves using the {asample} wildcard twice and then constructing the output directory name
> > to give to fastqc in the `-o` option.
> >
> > ~~~
> > rule fastqc:
> >   output:
> >     html = "{indir}.fastqc.{asample}/{asample}_fastqc.html",
> >     zip  = "{indir}.fastqc.{asample}/{asample}_fastqc.zip",
> >   input: "{indir}/{asample}.fq"
> >   shell:
> >     "fastqc -o {wildcards.indir}.fastqc.{wildcards.asample} {input}"
> > ~~~
> >
> {: .solution}
{: .challenge}

Our next option is to not worry about the individual files at all. Snakemake allows any output of a rule to be a
directory; you just have to specify this with the `directory()` function.
In this case, Snakemake doesn't track the individual files in the directory at all, so it makes the rule simpler.

We'll amend the `fastqc` rule so that the output is a single directory like so:

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{asample}")
  input:  "{indir}/{asample}.fq"
  shell:
     "fastqc -o {output} {input}"
~~~

On running this, we get an error.

~~~
$ snakemake -j1 -p reads.fastqc2.ref_1_1
...
Specified output directory 'reads.fastqc2.ref_1_1' does not exist
...
~~~

FastQC requires that the output directory must exist. (Other programs might insist that the output directory
does *not* exist.) The error can be rectified by making the directory explicitly in the `shell` code.

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{asample}")
  input:  "{indir}/{asample}.fq"
  shell:
     "mkdir {output} ; fastqc -o {output} {input}"
~~~

This works because the `shell` part of the rule can contain a whole script, with multiple commands to be run. Above
we used a semicolon to split the commands. For putting multiple lines into a `shell` section we use a special quoting syntax.

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{asample}")
  input:  "{indir}/{asample}.fq"
  shell:
     r"""mkdir {output}
         fastqc -o {output} {input}
      """
~~~

The "triple quoting" syntax comes from Python. Not only does it allow multiple lines to be added within the quotes but it
also allows you to embed both single and double quotes into the shell commands. The `r` character before the quotes disables
interpretation of "backslash escapes" like "\n" and "\t". This is good, as you want the Bash shell, not Snakemake itself, to
interpret these special characters. So when adding more complex shell sections, *always format them like this.*

This rule is also fine, but because the individual files are not explicitly named as outputs we may have problems chaining
later rules. Also consider that some applications won't give you any control over the output filenames. The most
powerful solution is to use shell commands to move and/or rename the files to exactly the names you want.

> ## Exercise - fixing FastQC to use our own output file names
>
> Complete the rule below so that the output filenames are correctly produced. You will need to add extra commands to the
> `shell` part aside from multiqc. Do not alter the `output` or `input` parts of the rule.
>
> ~~~
> rule fastqc:
>     output:
>         html = "{indir}.{asample}.fastqc.html",
>         zip  = "{indir}.{asample}.fastqc.zip"
>     input:  "{indir}/{asample}.fq"
>     shell:
>        r"""???
>         """
> ~~~
>
> > ## Solution
> >
> > This is one solution, using `-o .` to tell FastQC to output the files in the current directory.
> > They are then renamed to match the declared outputs. Remember that, after Snakemake runs all the shell commands,
> > it checks to see that all output file really were created.
> >
> > ~~~
> > rule fastqc:
> >     output:
> >         html = "{indir}.{asample}.fastqc.html",
> >         zip  = "{indir}.{asample}.fastqc.zip"
> >     input:  "{indir}/{asample}.fq"
> >     shell:
> >        r"""fastqc -o . {input}
> >            mv {wildcards.asample}_fastqc.html {output.html}
> >            mv {wildcards.asample}_fastqc.zip  {output.zip}
> >         """
> > ~~~
> {: .solution}
{: .challenge}

## Adding a Salmon rule

We saw above that the output of a rule can be a directory and saw the `directory()` function which declares this.
If you remember the rule for `kallisto quant` you may be thinking that this could have been written with
the whole directory as the output, and you would be right.

~~~
# Existing rule for kallisto_quant
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    ...

# Alternative rule with directory() output
rule kallisto_quant:
    output:  directory("kallisto.{sample}")
    ...
~~~

> ## Directories as input
>
> Note that you only use the `directory()` syntax for outputs. The input to a rule may be a directory without the need for any special
> syntax.
{: .callout}

> ## Exercise
>
> An alternative to kallisto for transcript quantification is `salmon`. The procedure is virtually identical, having an indexing
> step and a quantification step. Note that in real usage you are advised to prepare and add decoy sequences to the index, but for the
> purposes of this tutorial we'll just keep things as simple as possible.
>
> Based upon the following commands:
>
> ~~~
> $ salmon index -t <transcriptome as fastq> -i <index name> -k 31
> $ salmon quant -i <index name> -l A -1 <fastq1> -2 <fastq2> --validateMappings -o <output path>
> ~~~
>
> Add a pair of rules to index and align the reads with Salmon. Note that:
>
> 1. Unlike Kallisto, the index produced by Salmon is a directory of files, not a single file
> 1. As per the note above, you only need the `directory()` marker on the outputs of rules
>
{: .challenge}


## Combining the outputs with MultiQC

MultiQC scans for analysis outputs in a given directory and all subdirectories, then makes a report. It knows about FastQC, Salmon
and Kallisto outputs so we should be able to compile a report on all these. To scan the current directory, simply run:

~~~
$ multiqc . -o multiqc_out
~~~

> ## Exercise - adding a MultiQC rule
>
> Earlier we made a basic summary-type rule called "all_counts". Now make a "multiqc" rule that gathers up all the FastQC, Salmon
> and Kallisto reports.
>
> Considerations:
>
> 1. Your rule is going to have several named inputs, and these inputs will be lists of files generated with `expand()` functions.
> 1. Ensure that both "kallisto quant" and "salmon quant" are run on all 9 samples, that is all three repeats of all three conditions.
> 1. Since multiqc scans for input files, the input names don't have to be explicitly mentioned in the `shell` part.
>
> > ## Solution
> >
> > TODO
> >
> {: .solution}
>
> ### Bonus exercise
>
> Because `multiqc` scans for suitable input rather than taking an explicit list of files, there is a risk that it picks up
> unwanted reports if you leave old files in the directory. To make the rule fully explicit, one idea might be to make a
> temporary directory and symlink all the files into it, and then tell multiqc to look in there. Amend the rule so it does this.
>
> Tip - when making links, use `ln -sr <src> <target>` to avoid link relativity issues.
>
{: .challenge}

> ## Exercise - fixing Kallisto
>
> You may notice that MultiQC is not capturing Kallisto output when making the reports. The reason for this is given in the
> [MultiQC manual here](https://multiqc.info/docs/#kallisto):
>
> ~~~
> Note - MultiQC parses the standard out from Kallisto, not any of its output files (abundance.h5, abundance.tsv, and run_info.json).
> As such, you must capture the Kallisto stdout to a file when running to use the MultiQC module.
> ~~~
>
> Fix the Snakefile so that Kallisto standard output is redirected to a file and can be collected by MultiQC. (Hint - MultiQC does not
> mind what you call the file so choose your own sensible name).
>
> > ## Solution
> >
> > TODO
> {: .solution}
{: .challenge}

[fig-workflow]: ../fig/overview.svg

{% include links.md %}

