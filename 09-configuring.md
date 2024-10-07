---
title: Configuring workflows
teaching: 30
exercises: 20
---

::::::::::::::::::::::::::::::::::::::: objectives

- Add parameters to rules
- Use configuration files and command line options to set the parameters

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I separate my rules from my configuration?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the final Snakefile from episodes 1 to 8](files/ep08.Snakefile) you
may use to start this episode.*

## Adding parameters (params) to rules

So far, we've written rules with `input`, `output` and `shell` parts. Another useful section you
can add to a rule is `params`.

Consider the "trimreads" rule we defined earlier in the course.

```source
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
```

Can you remember what the `-t 20` and `-l 100` parameters do without referring back to the manual?
Probably not!
Adding comments in the Snakefile would certainly help, but we can also make important settings
into parameters.

```source
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    params:
        qual_threshold = "20",
        min_length     = "100",
    shell:
        "fastq_quality_trimmer -t {params.qual_threshold} -l {params.min_length} -o {output} <{input}"
```

Now it is a little clearer what these numbers mean. Use of parameters doesn't give you extra
functionality but it is good practise put settings like these into parameters as it makes the whole
rule more readable.

:::::::::::::::::::::::::::::::::::::::  challenge

## Adding a parameter to the `salmon_index` rule

Modify the existing salmon\_index rule so that the `-k` setting (k-mer length) is a parameter.

Change the length to 29 and re-build the index with the amended rule.

:::::::::::::::  solution

## Solution

```source
rule salmon_index:
    output:
        idx = directory("{strain}.salmon_index")
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    params:
        kmer_len = "29"
    shell:
        "salmon index -t {input.fasta} -i {output.idx} -k {params.kmer_len}"
```

```bash
$ snakemake -j1 -p -f Saccharomyces_cerevisiae.R64-1-1.salmon_index
```

Notes:

- You can choose a different parameter name, but it must be a valid identifier - no spaces or
  hyphens.
- Changing the parameters doesn't automatically trigger Snakemake to re-run the rule so you
  need to use `-f` (or `-R` or `-F`) to force the job to be re-run (but, as mentioned in
  [episode 5](05-the_dag.md),
  this behaviour is changed in recent Snakemake versions).

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Making Snakefiles configurable

It's good practise to break out parameters that you intend to change into a separate
file. That way you can re-run the pipeline on new input data, or with alternative settings, but
you don't need to edit the Snakefile itself.

We'll save the following lines into a file named *config.yaml*.

```source
salmon_kmer_len: "31"
trimreads_qual_threshold: "20"
trimreads_min_length: "100"
```

This file is in YAML format. This format allows you to capture complex data structures but we'll
just use it to store some name+value pairs. We can then reference these values within the
Snakefile via the **config** object.

```source
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    params:
        qual_threshold = config["trimreads_qual_threshold"],
        min_length     = config.get("trimreads_min_length", "100"),
    shell:
        "fastq_quality_trimmer -t {params.qual_threshold} -l {params.min_length} -o {output} <{input}"
```

In the above example, the **trimreads\_qual\_threshold** value must be supplied in the config, but
the **trimreads\_min\_length** can be omitted, and then the default of "100" will be used.

If you are a Python programmer you'll recognise the syntax here. If not, then just take note that
the first form uses *square brackets* and the other uses `.get(...)` with *regular brackets*. Both
the config entry name and the default value should be in quotes.

:::::::::::::::::::::::::::::::::::::::::  callout

## Using config settings with params

Strictly speaking, you don't have to use *config* in conjunction with *params* like this, but it's
normally a good idea to do so.

::::::::::::::::::::::::::::::::::::::::::::::::::

The final step is to tell Snakemake about your config file, by referencing it on the command line:

```bash
$ snakemake --configfile=config.yaml ...
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Making the parameter of a rule configurable

Fix the `salmon_index` rule to use `salmon_kmer_len` as in the config file sample above.
Use a default of "29" if no config setting is supplied.

Run Snakemake in *dry run* mode (`-n`) to check that this is working as expected.

:::::::::::::::  solution

## Solution

Rule is as before, aside from:

```source
params:
    kmer_len = config.get("salmon_kmer_len", "29")
```

If you run Snakemake with the `-n` and `-p` flags and referencing the config file, you should
see that the command being printed has the expected value of *31*.

```bash
$ snakemake -n -p -f --configfile=config.yaml Saccharomyces_cerevisiae.R64-1-1.salmon_index
```

*Note that if you try to run Snakemake with no config
file you will now get a **KeyError** regarding **trimreads\_qual\_threshold**. Even though you
are not using the **trimreads** rule, Snakemake needs a setting for all mandatory parameters.*

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

Before proceeding, we'll tweak the Snakefile in a couple of ways:

1. Set a default `configfile` option so we don't need to type it on every command line.
2. Get Snakemake to print out the config whenever it runs.

Add the following lines right at the top of the Snakefile.

```source
configfile: "config.yaml"
print("Config is: ", config)
```

Finally, as well as the `--configfile` option to Snakemake there is the `--config` option which
sets individual configuration parameters.

```bash
$ snakemake --configfile=config.yaml --config salmon_kmer_len=23 -p -nf Saccharomyces_cerevisiae.R64-1-1.salmon_index/
```

This is all getting quite complex, so in summary:

- Snakemake loads the `--configfile` supplied on the command line, or else defaults to the one
  named in the Snakefile, or else runs with no config file.
- Individual `--config` items on the command line always take precedence over settings in the
  config file.
- You can set multiple `--config` values on the command line and the list ends when there is
  another parameter, in this case `-p`.
- Use the `config.get("item_name", "default_val")` syntax to supply a default value which takes
  lowest precedence.
- Use `config["item_name"]` syntax to have a mandatory configuration option.

:::::::::::::::::::::::::::::::::::::::  challenge

## Making the conditions and replicates into configurable lists

Modify the *Snakefile* and *config.yaml* so that you are setting the *CONDITIONS* and
*REPLICATES* in the config file.
Lists in YAML use the same syntax as Python, with square brackets and commas, so you can copy the
lists you already have. Note that you're not expected to modify any rules here.

Re-run the workflow to make a report on *just replicates 2 and 3*. Check the MultiQC report to
see that it really does have just these replicates in there.

:::::::::::::::  solution

## Solution

In *config.yaml* add the lines:

```source
conditions: ["etoh60", "temp33", "ref"]
replicates: ["1", "2", "3"]
```

In the *Snakefile* we can reference the *config* while setting the global variables. There are
no *params* to add because these settings are altering the selection of jobs to be added to the
DAG, rather than just the *shell* commands.

```source
CONDITIONS = config["conditions"]
REPLICATES = config["replicates"]
```

And for the final part we can either edit the *config.yaml* or override on the command line:

```bash
$ snakemake -j1 -f --config replicates=["2","3"] -p multiqc_out
```

Note that we need to re-run the final report, but only this, so only `-f` is necessary. If you
find that replicate 1 is still in you report, make sure you are using the final version of the
*multiqc* rule from the previous episode, that symlinks the inputs into a
*multiqc\_in* directory.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is a Snakefile](files/ep09.Snakefile) incorporating the changes made in
this episode.*

*To run it, you need to save out [this file](files/ep09/config.yaml) to the same directory.*



:::::::::::::::::::::::::::::::::::::::: keypoints

- Break out significant options into rule parameters
- Use a YAML config file to separate your configuration from your workflow logic
- Decide if different config items should be mandatory or else have a default
- Reference the config file in your Snakefile or else on the command line with `--configfile`
- Override or add config values using `--config name1=value1 name2=value2` and end the list with a new parameter, eg. `-p`

::::::::::::::::::::::::::::::::::::::::::::::::::


