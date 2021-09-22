---
title: "Conda Integration"
teaching: 20
exercises: 10
questions:
- "How do I get Snakemake to manage software dependencies?"
objectives:
- "Understand how Snakemake works with Conda and Bioconda"
- "Add a new package to the workflow with Conda integration"
keypoints:
- "Add key points"
---

## Conda and Snakemake

### A recap on Conda

This episode assumes you have some familiarity with the basics of Conda (https://conda.io) and Bioconda (https://bioconda.github.io/).

To recap some key terms:

 * **Conda** is a system for installing software **packages** into self-contained directories called **environments**
 * **Conda** finds new packages in on-line repositories which are known as **channels**
 * The **Bioconda project** maintains a channel with many open source bioinformatics packages available
 * Each environment you create will have multiple packages, but if you want to install two versions of the same package, or two packages which
   otherwise conflict, you need to put them in **separate environments**
 * You can switch between environments using the **conda activate** command. This only affects the current shell session
 * An environment may be **exported**, which simply means asking Conda to list all the packages in that environment, as well as the channels
   where the packages are available, in a YAML format

Some key conda commands:

 * Make a new environment [with initial packages]
    * `$ conda create -n my-environment [python]`
 * List the environments currently known to conda
    * `$ conda env list`
 * Switch to using an environment in the current shell session
    * `$ conda activate my-environment`
 * Add a new default channel to all environments [or just the current one]
    * `$ conda config [--env] --append channels bioconda`
 * Install a single package into the current environment
    * `$ conda install snakemake`
 * Save a list of all packages in one environment in YAML format, then recreate an identical environment from that file.
    * `$ conda env export -n my-environment > my-environment.yaml`
    * `$ conda env create -p ~/cloned-environment -f my-environment.yaml`

> ## Exercises
>
> 1. Find out what version of *fastx_toolkit* is installed in the current conda environment.
>
> 1. Create a new conda env named *new-env* and install the *cutadapt* package from Bioconda into it.
>
> > ## Solution
> >
> > 1. There are several ways of doing this, but one is:
> >
> > ~~~
> > $ conda env export | grep fastx
> >  - fastx_toolkit=0.0.14=0
> > ~~~
> >
> > 2. Again there are different ways to achieve this but you can do it in a single command:
> >
> > ~~~
> > $ conda create -n new-env --channel bioconda --channel conda-forge cutadapt
> > ~~~
> >
> > This makes the environment and installs the latest cutadapt, searching in the *bioconda* and *conda-forge*
> > channels. As documented on the [Bioconda website](https://bioconda.github.io/user/install.html#set-up-channels),
> > some of the packages in Bioconda rely on other packages in Conda Forge, and this includes cutadapt, so you do need both.
> >
> > Depending on your system configuration, you may already have bioconda and conda-forge set as default
> > channels (use `conda config --show channels` to check), but in any case it doesn't hurt to include them explicitly.
> >
> {: .solution}
{: .challenge}

### Using conda with Snakemake

Up to now, the `shell` commands you have been running with Snakemake are whatever is installed on the local machine. The commands may
have been installed with conda, or with the system package manager, or installed manually. If you move your Snakefile to another machine,
or share it with a colleague, you will need to make sure all the dependencies are installed. If the versions of the packages are not the
same on two systems, you may discover that the workflow breaks, or produces different results.

If you are familiar with Conda, you may think to install the dependencies for your workflow into a conda environment, then `conda env export`
that into a YAML file, which you can use to quickly set up the same environment on another machine.

~~~
$ conda env export -n new-env > new-env.yaml
~~~

Snakemake takes this one step further by integrating directly with Conda. This brings some nice features...

Firstly it allows you to specify different environments to use for different rules. This is really useful if different rules need
mutually incompatible packages.

To configure this, add a `conda` declaration to a rule:

~~~
rule a_conda_rule:
    conda:  "new-env.yaml"
    shell:
        "which cutadapt"
~~~

Note that the declaration refers to the exported YAML file, not the existing environment. Snakemake will create the environment for you.
We can actually run the above rule, even though it doens't produce any data. We need to add the `--use-conda` flag when running
Snakemake.

~~~
$ snakemake -j1 --use-conda
Building DAG of jobs...
Creating conda environment new-env.yaml...
Downloading and installing remote packages.
Environment for new-env.yaml created (location: .snakemake/conda/d7df5e24)
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	a_conda_rule
	1
Select jobs to execute...

[Tue Sep 21 16:57:59 2021]
rule a_conda_rule:
    jobid: 0

Activating conda environment: /home/zenmaster/carpentries/nextflow_rnaseq_training_dataset/.snakemake/conda/d7df5e24
/home/zenmaster/carpentries/nextflow_rnaseq_training_dataset/.snakemake/conda/d7df5e24/bin/cutadapt
[Tue Sep 21 16:58:00 2021]
Finished job 0.
1 of 1 steps (100%) done
Complete log: /home/zenmaster/carpentries/nextflow_rnaseq_training_dataset/.snakemake/log/2021-09-21T165713.679409.snakemake.log
~~~

This takes some time, because Snakemake is (as the log says) creating the new environment. But if we run it a second time it's much
quicker. The environment that Snakemake made is stored under `./.snakemake/conda/d7df5e24`, and the last part of this name is a hash
of the contents of the `new-env.yaml` file. As long as the file is unchanged, Snakemake will re-use the same environment, but if the
YAML file is edited at all Snakemake will see that the hash has changed and make a new environment.

We'll do something useful with cutadapt in the next episode.

> ## Challenge
>
> Going back to our RNA-Seq workflow, imagine we want to try running the analysis with an older version of Salmon,
> to see if the results are different.
> We don't want to mess with the version of Salmon that's currently installed, or change any other parts of the
> workflow other than *salmon_index* and *salmon_quant*. Define a new Conda environment that contains `salmon=1.2.1`
> and then use this for the two rules by adding appropriate `conda:` directives. Then run your amended workflow.
>
> > ## Solution
> >
> > The first thing that we need is an appropriate YAML file. It's possible to write one from scratch, but I'll
> > follow the process shown above for cutadapt, that is to *conda create* an environment and then *export* it.
> >
> > ~~~
> > $ conda create -n salmon-1.2.1 --channel bioconda --channel conda-forge salmon=1.2.1
> > $ conda env export -n salmon-1.2.1 > salmon-1.2.1.yaml
> > ~~~
> >
> > Having done this, we can delete the environment as we just need the exported YAML file.
> >
> > ~~~
> > $ conda env remove -n salmon-1.2.1
> > ~~~
> >
> > Next add the *conda* declarations to both the *salmon_index* and *salmon_quant* rules.
> >
> > ~~~
> > conda: "salmon-1.2.1.yaml"
> > ~~~
> >
> > And when running the workflow, we need to use the `--use-conda` flag as well as telling Snakemake that these two rules
> > have changed and must be re-run.
> >
> > ~~~
> > $ snakemake -j1 -Rsalmon_index -Rsalmon_quant --use-conda ...
> > ~~~
> >
> {: .solution}
{: .challenge}


{% include links.md %}

