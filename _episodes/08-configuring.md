---
title: "Configuration files and parameters"
teaching: 20
exercises: 20
questions:
- "How do I separate my rules from my configuration?"
objectives:
- "Add parameters to rules"
- "Use configuration files and command line params"
keypoints:
- "Add key points"
---

## Adding parameters (params) to rules

So far, we've written rules with `input`, `output` and `shell` parts. Another useful section you can add to
a rule is `params`.

Consider the "trimreads" rule we defined earlier in the course.

~~~
rule trimreads:
  output: "trimmed/{asample}.fq"
  input:  "reads/{asample}.fq"
  shell:
    "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
~~~

Can you remember what the `-t 20` and `-l 100` parameters do without referring back to the manual? Probably not!
Adding comments in the Snakefile might be a good idea, but we can also make important settings into parameters.

~~~
rule trimreads:
  output: "trimmed/{asample}.fq"
  input:  "reads/{asample}.fq"
  params:
    qual_threshold = "20",
    min_length     = "100",
  shell:
    "fastq_quality_trimmer -t {params.qual_threshold} -l {params.min_length} -o {output} <{input}"
~~~

Now it is a little clearer what these numbers mean. Use of parameters does not give you extra functionality but it is
good practise to make important settings into parameters as it makes the rule more readable.

> ## Exercise
>
> Modify the existing salmon_index rule so that the `-k` setting (k-mer length) is a parameter.
>
> Change the length to 33 and re-build the index with the amended rule.
>
> > ## Solution
> >
> > ~~~
> > rule salmon_index:
> >     output:
> >         idx = directory("{strain}.salmon_index")
> >     input:
> >         fasta = "transcriptome/{strain}.cdna.all.fa.gz"
> >     params:
> >         kmer_len = "33"
> >     shell:
> >         "salmon index -t {input.transcriptome} -i {output.index} -k {params.kmer_len}"
> > ~~~
> >
> > ~~~
> > snakemake -j1 -p -f Saccharomyces_cerevisiae.R64-1-1.salmon_index
> > ~~~
> >
> > Notes:
> >
> > * You can choose a different parameter name, but it must be a valid identifier, no spaces or hyphens.
> > * Changing the parameters does automatically trigger Snakemake to re-run the rule (remember it only looks
> >   at file modification times) so you need to use `-f` (or `-R` or `-F`) to force the job to be re-run.
> >
> {: .solution}
{: .challenge}

## Making Snakefiles configurable

In general, it's good practise to break out parameters that you intend to change into a separate file. That
way you can re-run the pipeline on new input data, or with alternative settings, but you don't
need to edit the Snakefile itself.

We'll save the following lines into a file named "config.yaml".

~~~
salmon_kmer_len: "33"
trimreads_qual_threshold: "20"
trimreads_min_length: "100"
~~~

This file is in YAML format. This format allows you to capture complex data structures but we'll just use it to
store some name/value pairs. We can then reference these values within the Snakefile.

~~~
rule trimreads:
  output: "trimmed/{asample}.fq"
  input:  "reads/{asample}.fq"
  params:
    qual_threshold = config["trimreads_qual_threshold"],
    min_length     = config.get("trimreads_min_length", "100"),
  shell:
    "fastq_quality_trimmer -t {params.qual_threshold} -l {params.min_length} -o {output} <{input}"
~~~

You don't have to use the config values in conjunction with params like this, but it's often a good idea to do so.
In the above example, the **trimreads_qual_threshold** value must be supplied in the config, but the
**trimreads_min_length** can be omitted and then the default of "100" will be used.

If you are a Python programmer you'll recognise the syntax here. If not, then just take note that one form uses square
brackets and the other uses `.get` and regular brackets.

The final step is to tell Snakemake to load your config file. You can do this on the command line:

~~~
$ snakemake --configfile config.yaml ...
~~~

> ## Exercise
>
> Fix the `salmon_index` rule to use `salmon_kmer_len` as in the config file sample above. Use a default of "31" if
> no config setting is supplied.
>
> > ## Solution
> >
> > TODO
> >
> {: .solution}
{: .challenge}

Before proceeding, we'll tweak the Snakefile in a couple of ways:

1. Set a default `configfile` option.
1. Get Snakemake to print out the config whenever it runs.

Add the following lines right at the top of the Snakefile.

~~~
configfile: "config.yaml"
print("Config is: ", config)
~~~

Finally, as well as the `--configfile` option to Snakemake there is the `--config` option which sets individual
configuration parameters.

~~~
$ snakemake --configfile config.yaml --config salmon_kmer_len=23
~~~

This is all getting quite complex, so in summary:

* Snakemake loads the `--configfile` supplied on the command line, or else defaults to the one named in the Snakefile, or else
  runs with no config.
* Individual `--config` items always take precedence over settings in the config file.
* Use the `config.get("item_name", "default_val")` syntax to supply a default value which takes lowest precedence.

> ## Exercise
>
> Modify the Snakefile and config.yaml so that you are setting the CONDITIONS and REPLICATES in the config file.
> Note that setting lists works just the same way as single values.
>
> Re-run the workflow to make a report on just replicates 2 and 3. Check the MultiQC report to see that it really
> does have just these replicates in there.
>
{: .challenge}
