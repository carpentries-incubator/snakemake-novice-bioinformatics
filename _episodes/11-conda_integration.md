---
title: "Conda Integration"
teaching: 20
exercises: 10
questions:
- "How do I get Snakemake to manage software dependencies?"
objectives:
- "Understand how Snakemake works with conda and Bioconda"
- "Add a new package to the workflow with conda integration"
keypoints:
- "Add key points"
---

## Conda and Snakemake

### A recap on Conda

This episode assumes you have some familiarity with the basics of Conda (https://conda.io) and Bioconda (https://bioconda.github.io/).

To recap some key terms:

 * **Conda** is a system for installing software packages into self-contained directories called **environments**
 * **Conda** finds new packages in on-line repositories which are known as **channels**
 * The **Bioconda project** maintains a channel with many open source bioinformatics packages available
 * Each environment you create may have multiple packages, but if you want to install two versions of the same package, or two packages which
   conflict, you need to put them in **separate environments**
 * You can switch between environments using the **conda activate** command. This only affects the current shell session
 * An environment may be **exported**, which simply means asking Conda to list all the packages in that environment, as well as the channels
   where the packages are available

Some key conda commands:

 * Make a new environment
    * `$ conda create -n my_environment`
 * List the environments currently known to conda
    * `$ conda env list`
 * Switch to using an environment in the current shell session
    * `$ conda activate my_environment`
 * Add a new channel to be available to all environments [or just the current one]
    * `$ conda config [--env] --append channels bioconda`
 * Install a single package into the current environment
    * `$ conda install snakemake`
 * Print a list of all packages in the current environment in YAML format
    * `$ conda env export -n my_environment`

> ## Exercises
>
> 1. Find out what version of *fastx_toolkit* is installed in the current conda environment.
>
> 1. Create a new conda env named *new_env_1* and install the *cutadapt* package from Bioconda into it.
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
> > $ conda create -n new_env --channel bioconda --channel conda-forge cutadapt
> > ~~~
> >
> > This makes the environment and installs the latest cutadapt, searching in the *bioconda* and *conda-forge*
> > channels. As documented on the [Bioconda website](https://bioconda.github.io/user/install.html#set-up-channels),
> > some of the packages in Bioconda
> > rely on other packages in Conda Forge, and this includes cutadapt, so you do need both.
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

Snakemake takes this one step further with the following features:

1) Snakemake allows you to specify different environments to use for different rules. This is really useful if different rules need incompatible packages.
2) Snakemake can set up the environment for you, based on the YAML file made by `conda env export`

Export the env you built above
Then mention simplified YAML file
Then show the required directives (in the Snakefile and on the command line) that correspond to the two things above

When cutadapt runs, add a `which cutadapt` so we see it is running in an environment. Yeah.

~~~
Simplified file here
~~~

Now let's make a very simple workflow that uses cutadapt on just one sequence. Yeah. That should do.

Make this an exercise that leads to the longer challenge. Does the longer challenge want an episode of it's own? Mmmmmaybe.
