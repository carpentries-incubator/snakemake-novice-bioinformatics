---
title: "Chaining rules"
teaching: 40
exercises: 30
questions:
- "How do I combine rules into a workflow?"
- "How do I make a rule with multiple inputs and outputs?"
objectives:
- "Use Snakemake to filter and then count the lines in a FASTQ file"
- "Add an RNA quantification step in the data analysis"
- "See how Snakemake deals with missing outputs"
keypoints:
- "Snakemake links rules by iteratively looking for rules that make missing inputs"
- "Rules may have several named inputs and/or outputs"
- "If a shell command does not yield an expected output then Snakemake will regard that job as failed"
---
*For reference, [this is the Snakefile](../code/ep02.Snakefile) you should have to start the episode.*

## A pipeline of multiple rules

We now have a "trimreads" rule and a "countreads" rule. Following the previous chapter, the contents
of the Snakefile should be:

~~~
# New generic read counter
rule countreads:
  output: "{sample}.fq.count"
  input:  "reads/{sample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# Trim any FASTQ reads for base quality
rule trimreads:
  output: "trimmed/{sample}.fq"
  input:  "reads/{sample}.fq"
  shell:
    "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
~~~
{: .language}

The problem is there is no good way to use these rules together, that is, to trim an input file and
then count the reads in the trimmed file. The *countreads* rule only takes input reads from the *reads*
directory, whereas the *trimreads* rule puts all results into the *trimmed* directory.

Chaining rules in Snakemake is a matter of choosing filename patterns that connect the rules. There's something of
an art to it - most times there are several options that will work. Consider the following alternative version
of the *countreads* rule:

~~~
# New even-more-generic read counter
rule countreads:
  output: "{indir}.{sample}.fq.count"
  input:  "{indir}/{sample}.fq"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"
~~~
{: .language}

Now, the rule no longer requires the input files to be in the "reads" directory. The directory name has been
replaced by the `{indir}` wildcard. We can request Snakemake to create a file following this output pattern:

~~~
$ snakemake -j1 -F -p trimmed.ref1_1.fq.count
$ cat trimmed.ref1_1.fq.count
14278
~~~
{: .language-bash}

Look at the logging messages that Snakemake prints in the terminal. What has happened here?

1. Snakemake looks for a rule to make `trimmed.ref1_1.fq.count`
1. It determines that "countreads" can make this if `indir=trimmed` and `sample=ref1_1`
1. It sees that the input needed is therefore `trimmed/ref1_1.fq`
<br/><br/>
1. Snakemake looks for a rule to make `trimmed/ref1_1.fq`
1. It determines that "trimreads" can make this if `sample=ref1_1`
1. It sees that the input needed is therefore `reads/ref1_1.fq`
<br/><br/>
1. Now Snakemake has reached an available input file, it runs both steps to get the final output

This, in a nutshell, is how we build workflows in Snakemake.

1. Define rules for all the processing steps
1. Choose `input` and `output` naming patterns that allow Snakemake to link the rules
1. Tell Snakemake to generate the final output files

If you are used to writing regular scripts this takes a little
getting used to. Rather than listing steps in order of execution, you are always **working
backwards** from the final desired result. The order of operations is determined by applying the
pattern matching rules to the filenames, not by the order of the rules in the Snakefile.

This logic of working backwards from the desired output is why we're putting the `output` lines first
in all out rules - to remind us that these are what Snakemake looks at first!

> ## Thinking about your own workflows
>
> Think about any data processing task you have done yourself, and write down three or four steps from that
> workflow.
>
> What were the inputs to, and outputs from, each step?
>
> How did the steps connect up, in terms of data going from one to the next? You may want to sketch
> this out and use arrows to indicate the linkages between the steps.
>
> > ## A sample answer based on brewing a mug of tea.
> >
> > TODO - draw this out too.
> > ~~~
> > "Boil Water" : input=cold water, output=hot water
> >
> > "Brew Tea" : input=[hot water, tea bag, mug], output=tea infusion
> >
> > "Add Milk And Sugar" : input=[tea infusion, milk, sugar], output=cuppa
> > ~~~
> {: .solution}
{: .challenge}

## Adding an alignment step to the pipeline

Let's add another rule to our Snakefile. The reads we have are from a yeast RNA-seq experiment so we
want to quantify transcript abundance using the **kallisto** aligner. The command to do so looks
like this:

~~~
$ kallisto quant -i index_file -o output_dir in_1.fastq in_2.fastq
~~~

This command has three input files:

1. The transcriptome index
1. The first of the paired FASTQ files
1. The second of the paired FASTQ files

And it produces a directory of output files. According to
[the Kallisto manual](https://pachterlab.github.io/kallisto/manual#quant) this directory will
have three output files in it:

1. abundances.h5
1. abundances.tsv
1. run_info.json

We'll not worry about what the contents of these files mean just now. We just know that we want to
run the `kallisto quant` command and have it make the output files.

Making a rule with multiple inputs and outputs like this works much like the previous rules.

~~~
rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundances.h5",
        tsv  = "kallisto.{sample}/abundances.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"
~~~
{: .language}

There are many things to note here:

1. The individual input and output files are given names using the `=` syntax
1. Each of these lines must end with a `,` (optional for the last one)
1. In the `shell` part, the input placeholders are now like `{input.name}`
1. We've chosen to always quantify the trimmed version of the reads
1. Because `kallisto quant` only takes the output directory name, we've used the placeholder
   `{wildcards.sample}` rather than `{output}` which would give the full file names
1. We don't actually have the `{input.index}` file yet. This will need to be created using the `kallisto index`
   command
1. If the number of input or output files had been variable, we'd need a slightly different approach. We'll
   come on to this in a later episode.

Even though the rule is not going to work without the index, we can still run it to check that Snakemake is
happy with the rule definition.

> ## Running Kallisto on all replicates
>
> If you know about the Kallisto software, you may be thinking about running Kallisto on all replicates of
> the sample at once. We'll look at this later in the course, but for now assume that Kallisto is run once
> for each pair of FASTQ files.
>
{: .callout}

> ## Running the kallisto_quant rule
>
> Given that the *index* input is missing, what would you expect Snakemake to do if the new rule was run now?
>
> Try it by telling Snakemake to run the new rule on the files `ref1_1.fq` and `ref1_2.fq`.
> Since the rule defines multiple outputs, asking for any one of the output files will be enough.
>
> > ## Solution
> >
> > ~~~
> > $ snakemake -j1 -F -p kallisto.ref1/abundances.h5
> > ~~~
> >
> > Resulting error is "Missing input files". Note that Snakemake does not try to run kallisto at all. It stops
> > while still making a work plan, because it sees that a required input file is missing.
> >
> {: .solution}
{: .challenge}


> ## Building the index
>
> Instruct Snakemake how to build the genome index as part of the pipeline by adding another rule. The command
> we need to run is:
>
> ~~~
> $ kallisto index -i index_file_to_make fasta_file_to_index
> ~~~
>
> The file to be indexed is `transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz`. As there is only
> one input to the rule you don't have to give it a name, but you may do so if you like.
>
> Make it so that the output printed by the program is captured to a file, and therefore your rule will have
> two separate outputs: the index file and the log file. Note that the program prints messages on STDERR, so
> you will need to use `>&` rather than `>` to capture the output.
>
> Also, once you get all this to run, you will still see an error after Snakemake runs the "kallisto_quant" step.
> We'll look at this error and how to remedy it next.
>
> > ## Solution
> >
> > The rule you write could look something like this, but there are many variations that will work just as well.
> > Since there is only one transcriptome in the project, you may feel that use of the `{strain}` wildcard
> > is overkill, but who's to say we might not want to use another in future?
> >
> > ~~~
> > rule kallisto_index:
> >     output:
> >         idx = "{strain}.kallisto_index",
> >         log = "{strain}.kallisto_log",
> >     input:
> >         fasta = "transcriptome/{strain}.cdna.all.fa.gz"
> >     shell:
> >         "kallisto index -i {output.idx} {input.fasta} >& {output.log}"
> > ~~~
> >
> {: .solution}
{: .challenge}

## Dealing with the "missing files" error

Once we get the new rules to run, you'll see something like this:

~~~
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
~~~

There's a lot to take in here. Some of the messages are very informative. Some less so.

1. Snakemake did actually run kallisto, as evidenced by the output from kallisto that we see
1. There is no obvious error message in the kallisto output
1. Snakemake complains some expected output files are missing: `kallisto.ref1/abundances.h5` and `kallisto.ref1/abundances.tsv`
1. The third output file `kallisto.ref1/run_info.json` was found but has now been removed by Snakemake
1. Snakemake suggest this might be due to "filesystem latency"

This last point is a red herring. "Filesystem latency" is not an issue here. We can investigate further by looking at
the `kallisto.ref1` subdirectory.

~~~
$ ls kallisto.ref1/
abundance.h5  abundance.tsv
~~~
{: .language-bash}

So the file names created by kallisto are not quite the same as we saw in the manual (note - the manual may have been fixed at the
point you are doing this course, but it was true back when the course was written!). Change the rule definition in the
`Snakefile` to use the correct names, then you should have everything working.

> ## Errors are normal
>
> Don't be disheartened if you see errors like the one above when first testing your new Snakemake pipelines. There is a lot that
> can go wrong when writing a new workflow, and you'll normally need several iterations to get things just right. One advantage of
> the Snakemake approach compared to regular scripts is that Snakemake fails fast when there is a problem, rather than ploughing on
> and potentially running junk calculations on partial or corrupted data. Another advantage is that when a step fails we can resume
> from where we left off, as we'll see in the next episode.
>
{: .callout}

Finally, check that the rules are still generic by processing the *temp33_1* sample:

~~~
$ snakemake -j1 -F -p kallisto.temp33_1/abundance.h5
~~~
{: .language-bash}

*For reference, [this is a Snakefile](../code/ep03.Snakefile) incorporating the changes made in this episode.*

{% include links.md %}

