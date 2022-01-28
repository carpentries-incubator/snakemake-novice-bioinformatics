---
title: "Handling awkward programs"
teaching: 40
exercises: 40
questions:
- "How do I handle real bioinformatics tools?"
- "How do I define a rule where the output is a directory?"
objectives:
- "Understand some different choices available when defining a rule"
- "Learn about the `directory()` function"
- "Add multiple commands to the shell section of a rule"
keypoints:
- "Different bioinformatics tools will have different quirks"
- "If programs limit your options for choosing input and output filenames, you have several ways to deal with this"
- "Use triple-quote syntax to make longer shell scripts with multiple commands"
- "Use these extra commands to set up inputs and rename outputs"
---
*For reference, [this is the Snakefile](../code/ep05.Snakefile) you should have to start the episode.*

We've now seen how to link rules in a pipeline and how to merge results at the final step. This is the basic
pattern for many analysis workflows. For simplicity, we used mostly basic shell commands, but we'll now replace
these with bioinformatics tools:

* **Salmon** is a alternative to Kallisto, using a different alignment algorithm
* **FastQC** calculates a variety of metrics on a FASTQ file and produces an HTML report and a ZIP file.
* **MultiQC** combines the reports from various tools, including FastQC, Kallisto, and Salmon, into a single HTML report.

## The full workflow we are constructing

![Our full QC workflow][fig-workflow]

Real programs like this can have quirks like:

* Not allowing you to specify *input* or *output* file locations
* Outputting a whole directory of files
* Creating temporary files
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

It's possible to supply multiple input files at once, but the resulting output is exactly the same as if processing the
files one at a time. Two files are produced for each FASTQ file, and these files appear in the same directory as the
input file.
The `fastqc` command does not let us choose the output filenames, but we can set the output directory name.

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
  output: "{indir}.{sample}.fq.count"
  input:  "{indir}/{sample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~

To do the same with FastQC, we have various options:

1. Work with the default file names produced by FastQC and leave the reports in the same directory with the FASTQ files.
1. Make the outputs in a new directory named, eg. "reads.fastqc.ref_1_1/" (similar to what we did with Kallisto).
1. Do this, but tell Snakemake that the *directory itself* is the output.
1. Force our preferred naming convention by renaming the FastQC output files within the rule.

We'll try all four, and see where this gets us.

### Option 1: Working with the default FastQC output files

> ## Exercise - adding a FastQC rule using the default output file names
>
> Fill in the **???** to make a working rule for FastQC where `indir` may be "reads" or "trimmed". Do not change the
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
> > Since the `shell` command is not to be changed, the output names will be dictated by fastqc.
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
> {: .solution}
{: .challenge}

This rule is fine, but maybe we don't want to put the reports in with the sequences. As a general principle, when
writing Snakemake rules, you want to be in charge of output file names. FastQC lets you specify the output
directory, so we can use that.

### Option 2: Using the default output filenames in a new directory

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
> > This involves using the {sample} wildcard twice and then constructing the output directory name
> > to give to fastqc in the `-o` option.
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

### Option 3: Using the directory() function

Our next option is to tell Snakemake not to worry about the individual files at all, and just say that the
output of the rule is a directory. This makes the rule definition simpler.

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
> You only use the `directory()` function for outputs. The input to a rule may be a directory without the need for any special
> syntax.
{: .callout}


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
  output: directory("{indir}.fastqc.{sample}")
  input:  "{indir}/{sample}.fq"
  shell:
     "mkdir {output} ; fastqc -o {output} {input}"
~~~

This works because the `shell` part of the rule can contain a whole script, with multiple commands to be run. Above
we used a semicolon to split the commands. For putting multiple lines into a `shell` section there is a special quoting syntax.

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
later rules. Also consider that some applications won't give you any control over the output filenames. The most
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
> > This is one solution, using `-o .` to tell FastQC to output the files in the current directory.
> > They are then renamed to match the declared outputs. Remember that, after Snakemake runs all the shell commands,
> > it checks to see that all output file really were created.
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
> {: .solution}
{: .challenge}

## Adding Salmon as an alternative to Kallisto

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

> ## Exercise - adding Salmon as an alternative to Kallisto
>
> An alternative to *kallisto* for transcript quantification is *salmon*. The procedure is virtually identical, having an indexing
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
> 1. Unlike Kallisto, the index produced by Salmon is a directory of files, not a single file - both
>    of these new rules will produce a directory as output.
> 1. As per the note above, you only need the `directory()` marker on the outputs of rules
>
> > ## Solution
> >
> > ~~~
> > rule salmon_quant:
> >     output: directory("salmon.{sample}")
> >     input:
> >         index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
> >         fq1   = "trimmed/{sample}_1.fq",
> >         fq2   = "trimmed/{sample}_2.fq",
> >     shell:
> >         "salmon quant -i {input.index} -l A -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output}"
> >
> > rule salmon_index:
> >     output:
> >         idx = directory("{strain}.salmon_index")
> >     input:
> >         fasta = "transcriptome/{strain}.cdna.all.fa.gz"
> >     shell:
> >         "salmon index -t {input.fasta} -i {output.idx} -k 31"
> > ~~~
> >
> > If you copied the *kallisto_index* rule and logged the output of *salmon_index* to a file this is fine.
> > Just make sure when copying and pasting that you change all the parts that need to change!
> >
> {: .solution}
{: .challenge}


## Combining the outputs with MultiQC

MultiQC scans for analysis outputs in a given directory and all subdirectories, then makes a report. It knows about FastQC, Salmon
and Kallisto outputs so we should be able to compile a report on all these. To scan the current directory, simply run:

~~~
$ multiqc . -o multiqc_out
~~~

> ## Exercise - adding a MultiQC rule
>
> Earlier we made a basic summary-type rule called *all_counts*. Now make a *multiqc* rule that gathers up all the FastQC, Salmon
> and Kallisto reports.
>
> Considerations:
>
> 1. Your rule is going to have several named inputs, and these inputs will be lists of files generated with `expand()` functions.
> 1. Ensure that both *kallisto_quant* and *salmon_quant* are run on all 9 samples, that is all three repeats of all three conditions.
> 1. Run FastQC on the untrimmed reads only, and note that MultiQC specifically uses the .zip files for input.
> 1. Since multiqc scans for input files, the input names don't have to be explicitly mentioned in the `shell` part.
>
> > ## Solution
> >
> > *Note: This answer assumes that the **kallisto_quant** rule has been modified with a directory output as above.*
> >
> > ~~~
> > rule multiqc:
> > output: directory('multiqc_out')
> > input:
> >     salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >     kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >     fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
> > shell:
> >     "multiqc . -o multiqc_out"
> > ~~~
> >
> {: .solution}
{: .challenge}

> ## Exercise - fixing Kallisto
>
> You may notice that MultiQC is not capturing Kallisto output when making the reports. The reason for this is given in the
> [MultiQC manual here](https://multiqc.info/docs/#kallisto):
>
> > *Note - MultiQC parses the standard out from Kallisto, not any of its output files (abundance.h5, abundance.tsv, and run_info.json).
> > As such, you must capture the Kallisto stdout to a file when running to use the MultiQC module.*
>
>
> Fix the Snakefile so that Kallisto terminal output is redirected to a file and can be collected by MultiQC.
>
>  * *Hint 1:* The manual above is not quite right - you need to capture both **stdout and stderr**, so using `>&` rather than `>`.
>  * *Hint 2:* MultiQC does not mind what you call the file, so choose your own sensible name.
>
> > ## Solution
> >
> > ~~~
> > # Kallisto quantification of one sample, with log capture.
> > rule kallisto_quant:
> >     output: directory("kallisto.{sample}")
> >     input:
> >         index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
> >         fq1   = "trimmed/{sample}_1.fq",
> >         fq2   = "trimmed/{sample}_2.fq",
> >     shell:
> >      r"""mkdir {output}
> >          kallisto quant -i {input.index} -o {output} {input.fq1} {input.fq2} 2> {output}/kallisto_quant.log
> >       """
> > ~~~
> >
> > There are several perfectly good ways of structuring this, so don't worry if your answer is different.
> >
> > A gotcha with the above version is that the output directory needs to be created before *kallisto quant* is run,
> > much like with FastQC.
> > Remember that Snakemake deletes any existing outputs, including outputs that are directories, before the job is
> > run, and while Kallisto will create the directory for you this is too late for the shell to make the log file
> > and without the `mkdir` command it will report "No such file or directory".
> >
> > Another option is to make the log file outside of the directory, or stick with declaring the individual files as
> > outputs, as when we first made the rule, in which case the directory will be made by Snakemake. It's possible to
> > declare both a directory and a file within that directory as separate outputs, but this is probably not the best
> > idea.
> >
> {: .solution}
{: .challenge}

> ## Bonus exercise - making the MultiQC rule more robust
>
> Because MultiQC scans for suitable input rather than taking an explicit list of files, there is a risk that it picks up
> unwanted reports if you leave old files sitting in the directory. To make the rule fully explicit, one idea is to make a
> temporary directory and symlink all the files into it, and then tell MultiQC to look in there. Amend the *multiqc* rule
> so it does this.
>
>  * *Hint 1:* When making links, use `ln -snr -t <target_dir> <src>` to avoid link relativity issues.
>  * *Hint 2:* If you feel that tweaking some other rules would make this easier, feel free to do that.
>
> > ## Solution
> >
> > This solution will work with the version of *kallisto_quant* in the solution above.
> >
> > ~~~
> > rule multiqc:
> >     output:
> >         mqc_out = directory('multiqc_out'),
> >         mqc_in  = directory('multiqc_in'),
> >     input:
> >         salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
> >     shell:
> >       r"""mkdir {output.mqc_in}
> >           ln -snr -t {output.mqc_in} {input}
> >           multiqc {output.mqc_in} -o {output.mqc_out}
> >        """
> > ~~~
> >
> {: .solution}
{: .challenge}

To see the actual MultiQC report, open the file *multiqc_out/multiqc_report.html* in a web browser. You can do this from the
command line:

~~~
$ xdg-open multiqc_out/multiqc_report.html
~~~
{: .language-bash}

The report has a few issues, but we'll not get distracted by the details of how to configure MultiQC to resolve them.

[fig-workflow]: ../fig/overview.svg

*For reference, [this is a Snakefile](../code/ep06.Snakefile) incorporating the changes made in this episode.
You may now proceed to any later episode using this workflow as a starting point.*

{% include links.md %}

