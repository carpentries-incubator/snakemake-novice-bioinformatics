---
title: Conda integration
teaching: 50
exercises: 20
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand how Snakemake works with Conda and Bioconda
- Modify a package in the workflow with Conda integration

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I install new packages with Conda?
- How do I get Snakemake to manage software dependencies?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the final Snakefile from episodes 1 to 7](code/ep07.Snakefile) you may
use to start this episode.*

## The basics of Conda

You may already have some familiarity with the basics of [Conda](https://conda.io) and
[Bioconda](https://bioconda.github.io/).

Some key terms:

- **Conda** is a system for installing software **packages** into self-contained directories
  called **environments**.
- **Conda** finds new packages in on-line repositories which are known as **channels**.
- The **Bioconda project** maintains a channel with many open source bioinformatics
  packages available.
- Old versions of these packages are available, as well as the latest builds.
- Any environment may have multiple packages, but to install two versions of the same package, or
  two packages which conflict, you need to put them in **separate environments**.
- You can switch between environments using the **conda activate** command.
- An environment may be **exported**, which simply means getting Conda to print all the packages
  in that environment, as well as the channels where those packages are available, in a
  YAML format.

Assuming you are running throught this course in the [recommended setup
](\(setup.md)({{ page.root }}{% link setup.md %}
), Conda is already set up on the systems we are using.

Some key conda commands:

- Make a new environment and give it a name
  - `$ conda create -n my-environment`
- List the environments currently known to conda
  - `$ conda env list`
- Switch to a named environment in the current shell session
  - `$ conda activate my-environment`
- Configure Conda to enable the Bioconda channel (the new settings will apply to all environments)
  - `$ conda config --add channels bioconda`
  - `$ conda config --add channels conda-forge`
- Install the latest version of a single package into the current environment
  - `$ conda install snakemake`
- Install a specific version of a single package into the current environment
  - `$ conda install snakemake=5.14.0`
- Save a list of all packages in one environment in YAML format, then recreate an identical
  environment from that information.
  - `$ conda env export -n my-environment > my-environment.yaml`
  - `$ conda env create -n cloned-environment -f my-environment.yaml`

:::::::::::::::::::::::::::::::::::::::::  callout

## Channel configuration and conda-forge

For advanced usage, there are many ways you might want to configure your Conda channels. It's
even possible to have specific channel settings for each environment. For our purposes, we just
want to have the *bioconda* channel working, and as noted on
[the Bioconda website](https://bioconda.github.io/user/install.html#set-up-channels) this
involves a dependency a second channel named *conda-forge*, which provides some supporting tools
and libraries.

To be sure all is well, check your channel settings:

```output
$ conda config --show channels
channels:
  - conda-forge
  - bioconda
  - defaults
```

If you don't see these exact channels in this order, try the `conda config ...` commands shown
above to fix the situation. Once this configuration is right you won't need to do anything else
regarding the channel configuration in this course.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Exercises

1. Find out what version of *fastx\_toolkit* is installed in the current conda environment.

2. Create a new conda environment named *new-env* and install the *cutadapt* package from
  Bioconda into it.

:::::::::::::::  solution

## Solution

1. There are several ways of doing this, but one using the commands above is:

```bash
$ conda env export | grep fastx
 - fastx_toolkit=0.0.14=0
```

2. Conda should already be configured to install Bioconda packages (see the callout above) so
  we can do this:

```bash
$ conda create -n new-env
$ conda activate new-env
$ conda install cutadapt
```

It's also possible to install packages at the same time as creating the environment, though
this wasn't shown in the examples above.

```bash
$ conda create -n new-env cutadapt
```

You'll still need to *activate* the new environment in order to run *cutadapt*.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Using Conda with Snakemake

Up to now, the `shell` commands you have been running via Snakemake have called programs in your
default *PATH*. The commands may have been installed with conda, or with the system package
manager, or installed manually. If you move your Snakefile to another machine, or share it with
a colleague, you will need to make sure all the dependencies are installed. If the versions of the
packages are not the same on the two systems, you may discover that the workflow breaks, or
produces different results.

Once you are familiar with Conda, you may think to install the dependencies for your workflow into
a Conda environment, then `conda env export` that into a YAML file, which you can use to quickly
set up the same environment on any other machine.

```bash
$ conda env export -n new-env > new-env.yaml
```

Snakemake takes this one step further by integrating directly with Conda. This brings some
nice features...

Firstly it allows you to specify different environments to use for different rules. This is really
useful if different rules need mutually incompatible packages.

To configure this, add a `conda` declaration to a rule:

```source
rule a_conda_rule:
    conda:  "new-env.yaml"
    shell:
        "which cutadapt"
```

Note that the declaration refers to the exported YAML file, not any existing environment. Snakemake
will create the environment for you. We can now run the above rule, even though it doesn't
produce any outputs. We need to add the `--use-conda` option when running Snakemake.

```output
$ snakemake -j1 --use-conda a_conda_rule
Building DAG of jobs...
Creating conda environment new-env.yaml...
Downloading and installing remote packages.
Environment for new-env.yaml created (location: .snakemake/conda/d7df5e24)
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
    count   jobs
    1       a_conda_rule
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
```

This takes some time, because Snakemake is (as the log says) creating the new environment and
installing the packages. But if we run it a second time it's much quicker. This brings us to a
second feature of Snakemake+Conda integration: as long as the YAML file is unchanged, Snakemake
will re-use the same environment, but if the file is edited or replaced Snakemake will detect the
change and make a new environment. This happens because the environment that Snakemake made is
stored under `./.snakemake/conda/d7df5e24`, and the last part of this name is a hash of the
contents of the `new-env.yaml` file.

Also, it may seem very inefficient to create the environment twice, but Conda employs file
de-duplication so you don't actually end up with two copies of all the software.

We'll do something useful with *cutadapt* in the next episode.

:::::::::::::::::::::::::::::::::::::::  challenge

## Challenge

Going back to our RNA-Seq workflow, imagine we want to try running the analysis with an older
version of Salmon, to see if the results are different. We'll use **Salmon 1.2.1** which is
available in Bioconda.

This particular package happens to need a dependency, `tbb=2020.2`, which is not installed by
default but we can explicitly install it into the environment with Salmon. Without this Salmon
1\.2.1 will crash out with the message *error while loading shared libraries*.

We don't want to mess with the version of Salmon that's currently installed, or change any parts
of the workflow other than adding directives to *salmon\_index* and *salmon\_quant*. Define a new
Conda environment that includes the packages `salmon=1.2.1` and `tbb=2020.2` and then use this
for the two rules by adding appropriate `conda:` directives. Then run your amended workflow.

:::::::::::::::  solution

## Solution

The first thing that we need is an appropriate YAML file. It's possible to write one from
scratch in a text editor, but for this answer we'll follow the process shown above, that is to
*conda create* an environment, install the required packages and then *export* it.

```bash
$ conda create -n salmon-1.2.1
$ conda activate salmon-1.2.1
$ conda install salmon=1.2.1 tbb=2020.2
$ conda env export -n salmon-1.2.1 > salmon-1.2.1.yaml
```

As an aside, having done this, we could delete the environment as we just need the exported
YAML file.

```bash
$ conda deactivate
$ conda env remove -n salmon-1.2.1
```

Next add the same *conda* declaration to both the *salmon\_index* and *salmon\_quant* rules.

```source
conda: "salmon-1.2.1.yaml"
```

And when running the workflow, we need to give the `--use-conda` option as well as telling
Snakemake that these two rules have changed and must be re-run.

```bash
$ snakemake -j1 -Rsalmon_index -Rsalmon_quant --use-conda ...
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Note

If you look into a YAML file created with `conda env export` you will see that Conda lists every
single package dependency and the list is quite large. You may prefer to write your own YAML
file where you can be as precise as you like about which packages are needed and which version or
versions are acceptable. Conda will fill in the blanks to ensure all dependencies are met.

A file for an environment with cutadapt could be as simple as this.

```source
name: cutadapt-env
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - cutadapt=3.4
```

In this case, Conda will install version 3.4 of Cutadapt, but will meet the required
dependencies by installing the newest packages found in the channels, so the exact
environment will not be the same each time. This may make the workflow more compatible if,
for example, you try to switch from Linux to a Mac, but it may also cause problems if the newer
packages somehow don't work properly with Cutadapt 3.4, as happens with the `tbb` package above.
Sadly, while Conda is a great tool for installing and managing software it does have quirks and
shortcomings, and software setup continues to be a perennial headache for bioinformaticians.

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is a Snakefile](code/ep10.Snakefile) incorporating the changes made in
this episode.*



:::::::::::::::::::::::::::::::::::::::: keypoints

- Conda is a system for managing software packages in self-contained environments
- Snakemake rules may be associated with specific Conda environments
- When run with the `--use-conda` option, Snakemake will set these up for you
- Conda gives you fine control over software versions, without modifying globally-installed packages
- Workflows made this way are super-portable, because Conda handles installing the correct versions of everything

::::::::::::::::::::::::::::::::::::::::::::::::::


