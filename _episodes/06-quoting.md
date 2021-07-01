---
title: "Handling awkward programs - directory output and wrappers"
teaching: 0
exercises: 0
questions:
- "How do I handle real bioinformatics tools?"
- "How do I make my Snakefiles robust?"
objectives:
- "Handle commands that don't quite behave"
- "Understand issues regarding quoting special characters"
keypoints:
- "Add key points"
---

We've now seen how to link rules in a pipeline and how to merge results at the final step. For simplicity, we made
a rule that counts sequences and a rule that combines the results using a simple `cat` command. We'll now replace
these with the FastQC and MultiQC tools.

* FastQC calculates a variety of metrics on a FASTQ file and produces an HTML report and a ZIP file.
* MultiQC combines the reports from various tools, including FastQC, Kallisto, and Salmon, into a single HTML report.

## Adding a FastQC rule

We'll first run FastQC on a single input file and see what is produced.

~~~
$ fastqc reads/ref_1_1.fq
$ ls -ltr reads
...
-rw-r--r-- 1 zenmaster  users   464852 Jun  9 14:31 ref_1_1_fastqc.zip
-rw-r--r-- 1 zenmaster  users   654810 Jun  9 14:31 ref_1_1_fastqc.html
~~~

It's possible to supply multiple input files, but the output is exactly the same as if processing the files one at
a time. Two files are produced for each FASTQ file, and these files appear in the same directory as the input file.
The `fastqc` command does not let us specify the output filenames, but we can set the output directory name (-o option).

For the `countreads` rule we wrote earlier, we chose our preferred output file names first, then wrote the shell command
to put the results into that filename:

~~~
# Our existing countreads rule:
rule countreads:
  output: "{indir}.{asample}.fq.count"
  input:  "{indir}/{asample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~

We have three main options here:

1. Work with the file names produced by FastQC and leave the reports in the same directory with the FASTQ files.
1. Make each report in a new directory named, eg. "reads.fastqc.ref_1_1/" (similar to what we did with Kallisto).
1. Force our preferred naming convention by renaming the FastQC output files.

All are valid options. We'll try all three.

> ## Exercise - adding a FastQC rule by the first approach
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
writing Snakemake rules, you want to be in charge of where files are written. FastQC lets you specify the output
directory, so we can use that.

## Did we actually do directories as outputs yet? Not sure, once its all rearranged ##

> ## Exercise - adding a FastQC rule by the second approach
>
> Amend the `fastqc` rule so that the output is a single directory like so:
>
> ~~~
> rule fastqc:
>   output: directory("{indir}.fastqc.{asample}")
> ~~~
>
> Hint - run `fastqc --help` to see options for FastQC.
>
> ## Solution
>
> ~~~
> rule fastqc:
>   output: directory("{indir}.fastqc.{asample}")
>   input:  "{indir}/{asample}.fq"
>   shell:
>      "fastqc -o {output} {input}"
> ~~~

This rule is also fine, but because the individual files are not explicitly named as outputs we may have problems chaining
later rules. You should also consider that some applications won't give you any control over the output filenames. The most
powerful solution is to use shell commands to rename the files to exactly what you want.

The `shell` part of the rule can contain a whole script - multiple commands to be run. When putting multiple lines into a
`shell` section we use a special quoting syntax.

~~~
rule fastqc:
    output:
        html = "{indir}.{asample}.fastqc.html",
        zip  = "{indir}.{asample}.fastqc.zip"
    input:  "{indir}/{asample}.fq"
    shell:
       r"""fastqc -o . {input}
           mv {wildcards.asample}_fastqc.html {output.html}
           mv {wildcards.asample}_fastqc.zip  {output.zip}
        """
~~~

The "triple quoting" syntax comes from Python. Not only does it allow multiple lines to be added within the quotes but it
also allows you to embed both single and double quotes into the shell commands. The `r` character before the quotes disables
interpretation of "backslash escapes" like "\n" and "\t". This is good, as you want the BAsh shell, not Snakemake itself, to
interpret these special characters. So when adding more complex shell sections, always format them like this.

The actual shell commands tell FastQC to output the files in the current directory, but then rename them to match the declared
outputs. Remember that if Snakemake runs a job and the output files do not get created by the shell command then it will report
an error.

## Directories as inputs and outputs

Maybe this goes in the next chapter?? Yes it does!

In the exercise below, we'll work with a program that produces a whole directory of files as output. We already saw this
for `kallisto quant` and in this case the directory contained three files, so we listed these as three outputs of the rule.

~~~
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    ...
~~~

There are two other valid approaches:

1. List just a subset of the output files as outputs of the rule. Snakemake does not care (in fact does not check) if the
   command produces other files too. So for the `kallisto_quant` rule we could have just said `output: "abundance.h5"` and
   the rule would work. The other outputs still get created but Snakemake does not consider them when linking rules.

1. Tell snakemake that the output of the rule is a directory. This way, Snakemake will not consider the files inside the
   directory at all. Do this by adding `directory(...)` around the output path.

~~~
rules kallisto_quant:
    output: directory("kallisto.{sample}")
    ...
~~~

Note that you only have to do this for outputs. The input to a rule may be a directory without the need for any special
syntax.

> ## Challenge
>
> Adapt the Snakefile to run 'kallisto quant' on all 9 samples, that is all three repeats of all three conditions. Rather than
> listing the inputs explicitly, use `glob_wildcards` to find them as demonstrated above.
>
> An alternative to kallisto for transcript quantification is `salmon`. The procedure is virtually identical, having an indexing
> step and a quantification step. Note that in real usage one is advised to prepare and add decoy sequences to the index but for the
> purposes of this tutorial we'll just keep things simple.
>
> Based upon the following commands:
>
> ~~~
> $ salmon index -t <transcriptome as fastq> -i <index name> -k 31
> $ salmon quant -i <index name> -l A -1 <fastq1> -2 <fastq2> --validateMappings -o <output path>
> ~~~
>

## Running MultiQC

MultiQC scans for analysis outputs in a given directory and all subdirectories, then makes a report. It knows about FastQC, Salmon
and Kallisto outputs so we should be able to compile a report on these. To scan the current directory, simply run:

"multiqc . -o multiqc_out"

> ## Exercise - adding a MultiQC rule
>
> Earlier we made a basic summary-type rule called "all_counts". Now make a "multiqc" rule that gathers up all the FastQC, Salmon
> and Kallisto reports.
>
> Considerations:
>
> 1. Your rule is going to have several named inputs, and these inputs will be lists of files.
> 1. The inputs and outputs don't have to be explicitly mentioned in the `shell`, but you may decide to do so.

> ## Exercise - fixing Kallisto
>
> You may notice that MultiQC is not capturing Kallisto output when making the reports. The reason for this is given in the
> [MultiQC manual here](https://multiqc.info/docs/#kallisto).
>
> ~~~
> Note - MultiQC parses the standard out from Kallisto, not any of its output files (abundance.h5, abundance.tsv, and run_info.json).
> As such, you must capture the Kallisto stdout to a file when running to use the MultiQC module.
> ~~~
>
> Fix the Snakefile so that Kallisto standard output is saved to a file and can be collected by MultiQC. (Hint - MultiQC does not
> mind what you call the file so choose a sensible name).

###################
# Probably an episode split here. Parameters and quoting seems a good rule. Then I can do MultiQC first.
# Note - the Carpentries way seems to be that we intro quoting issues as they become relevant to the
# job in hand. So maybe do it that way? I think intro parameters here and quoting follows, then
# configuration from that.

## Adding Parameters to rules

So far, we've written rules with `input`, `output` and `shell` parts. Another useful section you can add to
a rule is `parameters`.

# Yeah but doesn't this belong in chapter 8? Maybe I should finish the DAG first?

Note - we probably want to introduce parameters first. The 'k' example in my assembly was perfect for this
but I can work out something else. Maybe a subsampling number? Or a quality cutoff?
###################

Consider the following simple Bash shell command:

  $ echo How "cool" are    you?
  How cool are you?

The message is printed back, but not before the shell has interpreted the text as per the standard command-line
parsing rules

 * The quotes around '"cool"' have been removed
 * The extra spaces after 'are   ' have been ignored
 * More subtly, the last work will be interpreted as a glob pattern...

  $ touch youX youY
  $ echo How "cool" are you?
  How cool are youX youY

Note: if you have certain shell settings you may see a warning about the unmatched glob pattern.

Anyone who has done any amount of programming or shell scripting will have come across quoting issues in code. Bugs related
to quoting can be very troublesome. In Snakemake these are particularly complex because the "shell" part of each rule
undergoes three rounds of interpretation before it is actually run:

 1. The string is parsed according to Python quoting rules
 1. Placeholders (eg. {input} {output}) are then replaced
 1. The resulting string goes to the shell and is subject to all Bash parsing rules

In this episode we introduce the best practises for making your Snakefiles robust, and some simple rules to avoid most
mis-quoting complications.

Now have some exercises that illustrate the three points and introduce ways to deal with them:

r"""Strings like this""". Allows embedded newlines and literal \n \t or "quotes". Ie. Python does as little
as possible to the string. Good practise.

Double-escaping {{curlies}}. There is no good way to avoid this. Or is there? We can add them outside of the rule, maybe?
The AWK example is good here.

Shell parsing rules. Normally you want this because you want your pipes, redirects, shell variables, etc. to work.
If you choose filenames using only [0-9A-Z-a-z.\_] then you will always be "safe", but sometimes you have to deal with
awkwardly named files or parameters. If you want them unmolested. Use {foo:q}. Not '{foo}'!
Have an exercise that demonstrates why this is good.

How will I best demonstrate the power of {:q}?

letters = ['"', "'", "$ ", "\n\tx", "{:}"]
shell:
    "echo {letters:q}"

Hmmm. A little obscure. Well have a think.

{% include links.md %}

