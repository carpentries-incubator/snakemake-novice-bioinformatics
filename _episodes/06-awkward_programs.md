---
title: "Handling awkward programs"
teaching: 40
exercises: 40
questions:
- "How do I handle tools which don't let me specify output file names?"
- "How do I define a rule where the output is a directory?"
objectives:
- "Add the FastQC tool to the pipeline"
- "Understand some different choices available when defining a rule"
- "Learn about the `directory()` function"
- "Add multiple command lines to the `shell` section of a rule"
keypoints:
- "Different bioinformatics tools will have different quirks"
- "If programs limit your options for choosing input and output filenames, you have several ways to deal with this"
- "Use triple-quote syntax to make longer shell scripts with multiple commands"
---
*For reference, [this is the Snakefile](../code/ep05.Snakefile) you should have to start the episode.*

## Introducing FastQC

FastQC is a popular tool for scanning FASTQ files and producing a selection of quality plots. It's "fast" in
that it runs quickly, and also in that you can normally use the default options so there is no configuration
needed.

![Screenshot of a typical FastQC report specifically showing the per-base quality box-and-whisker plot][fig-fastqc]

The program can be run interactively or in batch mode, where it saves out results as an HTML file plus a ZIP file.
We'll obviously need to use the batch mode to include it as part of our workflow, and we'll see shortly that
this presents a minor problem.

In general, real bioinformatics programs like FastQC can have quirks like:

* Not allowing you to specify *input* or *output* file names
* Outputting a whole directory of files
* Creating temporary files in various places
* Not supporting robust error handling

With a little care we can handle all these problems while still keeping our workflow tidy.

## Adding a FastQC rule

We'll first run FastQC on a single input file and see what is produced.

~~~
$ fastqc reads/ref_1_1.fq
$ ls -ltr reads
...
-rw-r--r-- 1 zenmaster  users   464852 Jun  9 14:31 ref_1_1_fastqc.zip
-rw-r--r-- 1 zenmaster  users   654810 Jun  9 14:31 ref_1_1_fastqc.html
~~~

It's possible to supply multiple input files, but the resulting output is exactly the same as if processing the
files one at a time. Two files are produced for each FASTQ file, and these files appear in the same directory as the
input file.
The `fastqc` command does not let us choose the output filenames, but we can set an alternative output directory.

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

For the `countreads` rule we wrote earlier (see episode 3), we chose our preferred output file name pattern so as to
allow us to effectively link rules. This gives us a rule that can count both *trimmed* and *untrimmed* reads.

~~~
# Our existing countreads rule, for reference...
rule countreads:
  output: "{indir}.{sample}.fq.count"
  input:  "{indir}/{sample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~

To do the same with FastQC, to report on both the trimmed and untrimmed reads, we have various options:

1. Work with the default file names produced by FastQC and leave the reports in the same directory with the FASTQ files.
1. Make the outputs in a new directory named, eg. "reads.fastqc.ref_1_1/" (similar to what we did with Kallisto).
1. Do this, but tell Snakemake that the *directory itself* is the output.
1. Force our preferred naming convention by renaming the FastQC output files within the rule.

We'll try all four, and see where this gets us.

### Option 1: Working with the default FastQC output files

> ## Exercise - adding a FastQC rule using the default output file names
>
> Fill in the `???` to make a working rule for FastQC where `indir` may be "reads" or "trimmed". Do not change the
> shell command or input pattern at all. Remember FastQC always makes two output files, so add two named outputs.
>
> ~~~
> rule fastqc:
>   output:
>       ???
>   input:  "{indir}/{sample}.fq"
>   shell:
>       "fastqc {input}"
> ~~~
>
> > ## Solution
> >
> > Since the `shell` command is not to be changed, the output names will be dictated by FastQC as we saw when
> > running the command directly in the terminal.
> >
> > ~~~
> > rule fastqc:
> >   output:
> >       html = "{indir}/{sample}_fastqc.html",
> >       zip  = "{indir}/{sample}_fastqc.zip",
> >   input:  "{indir}/{sample}.fq"
> >   shell:
> >       "fastqc {input}"
> > ~~~
> >
> > This rule contains wildcards, so in order to run it you specify one or more target output files:
> >
> > ~~~
> > $ snakemake -j1 -p reads/etoh60_1_1_fastqc.html reads/etoh60_1_2_fastqc.html
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

This rule is fine, but maybe we don't want to put the reports in with the sequences. As a general principle, when
writing Snakemake rules, we want to be in charge of output file names. FastQC lets us specify the output
directory, so we can use that...

### Option 2: Using the default output filenames in a new directory

> ## Exercise - a FastQC rule where the output files go into a new directory
>
> Modify the rule so that the output files go into a new directory. This will be very similar to the rule for
> `kallisto quant` added in episode 3.
>
> For example, when running on the file "trimmed/ref_1_1.fq" the outputs should be
> ~~~
> trimmed.fastqc.ref_1_1/ref_1_1_fastqc.html
> trimmed.fastqc.ref_1_1/ref_1_1_fastqc.zip
> ~~~
>
> > ## Solution
> >
> > This involves using the {sample} wildcard twice and then constructing the output directory name
> > to place in the `-o` option to fastqc.
> >
> > ~~~
> > rule fastqc:
> >   output:
> >     html = "{indir}.fastqc.{sample}/{sample}_fastqc.html",
> >     zip  = "{indir}.fastqc.{sample}/{sample}_fastqc.zip",
> >   input: "{indir}/{sample}.fq"
> >   shell:
> >     "fastqc -o {wildcards.indir}.fastqc.{wildcards.sample} {input}"
> > ~~~
> >
> {: .solution}
{: .challenge}

### Option 3: Using a `directory()` output

Our next option is to tell Snakemake not to worry about the individual files at all, and just say that the
output of the rule is a whole directory. This makes the rule definition simpler.

We'll amend the *fastqc* rule like so:

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{sample}")
  input:  "{indir}/{sample}.fq"
  shell:
     "fastqc -o {output} {input}"
~~~

> ## Note
>
> You only use the `directory()` declaration for outputs. Any input to a rule may be a directory without
> the need for any special syntax.
{: .callout}


On running this, we get an error.

~~~
$ snakemake -j1 -p reads.fastqc.ref_1_1
...
Specified output directory 'reads.fastqc.ref_1_1' does not exist
...
~~~

This error is being printed by FastQC. FastQC requires that the output directory must exist.
(Other programs might insist that the output directory does *not* already exist.)
The error can be rectified by making the directory explicitly in the `shell` code.

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{sample}")
  input:  "{indir}/{sample}.fq"
  shell:
     "mkdir {output} ; fastqc -o {output} {input}"
~~~

This works because the `shell` part of the rule can contain a whole script, with multiple commands to be run. Above
we used a semicolon to split the commands. For putting multiple lines into a `shell` section there is a special
quoting syntax.

~~~
rule fastqc:
  output: directory("{indir}.fastqc.{sample}")
  input:  "{indir}/{sample}.fq"
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
later rules. Also consider that some applications won't give you any control at all over the output filenames. The most
powerful solution is to use shell commands to move and/or rename the files to exactly the names you want.

### Option 4: Insisting on our own file names

> ## Exercise - fixing FastQC to use our own output file names
>
> Complete the rule below so that the output filenames are correctly produced. You will need to add extra commands to the
> `shell` part after running `fastqc`. Do not alter the `output` or `input` parts of the rule.
>
> ~~~
> rule fastqc:
>     output:
>         html = "{indir}.{sample}_fastqc.html",
>         zip  = "{indir}.{sample}_fastqc.zip"
>     input:  "{indir}/{sample}.fq"
>     shell:
>        r"""???
>         """
> ~~~
>
> > ## Solution
> >
> > This is one solution, using `-o .` to tell FastQC to produce the files in the current directory,
> > then explicitly renaming them to match the declared rule outputs. Remember that, after Snakemake runs
> > all the commands in the `shell` block, it checks to see that all the expected `output` files are
> > present, but within the block you can run any commands and make any files you like.
> >
> > ~~~
> > rule fastqc:
> >     output:
> >         html = "{indir}.{sample}_fastqc.html",
> >         zip  = "{indir}.{sample}_fastqc.zip"
> >     input:  "{indir}/{sample}.fq"
> >     shell:
> >        r"""fastqc -o . {input}
> >            mv {wildcards.sample}_fastqc.html {output.html}
> >            mv {wildcards.sample}_fastqc.zip  {output.zip}
> >         """
> > ~~~
> >
> {: .solution}
{: .challenge}

> ## Note
>
> There is actually a problem with the above solution which only starts to matter when we allow Snakemake
> to run multiple jobs in parallel. Right now we are always using `-j1`, but if we used, eg. `-j2`,
> then potentially Snakemake may try to make "reads.ref_1_1_fastqc.html" and
> "trimmed.ref_1_1_fastqc.html" in parallel, and both instances would be trying to write to the same
> temporary files at the same time. Snakemake has an elegant solution to this, in the form of `shadow`
> rules, but we're getting ahead of ourselves. For now we're running one job at a time, and this will work.
>
{: .callout}

## Altering the Kallisto rule to have a `directory()` output

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

Make the change in your Snakefile now. In other workflows this might not be the right approach but in this case it
works fine and makes the Snakefile neater. It will also make sense in the next chapter where we add the remaining
rules and finish the workflow.

[fig-fastqc]: ../fig/fastqc.png

*For reference, [this is a Snakefile](../code/ep06.Snakefile) incorporating the changes made in this episode.*

{% include links.md %}
