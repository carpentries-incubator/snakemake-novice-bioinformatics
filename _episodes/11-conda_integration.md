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

This episode assumes you have some familiarity with the basics of Conda (https://conda.io) and Bioconda (https://bioconda.github.io/)

To recap some key terms:

 * conda is a system for installing software packages into self-contained directories called **environments**
 * conda finds new packages in on-line repositories which are known as **channels**
 * The **Bioconda project** maintains a  channel with many open source bioinformatics packages available
 * each environment you create may have multiple packages, but if you want to install two versions of the same package, or two packages which
   conflict, you need to put them in **separate environments**
 * you can switch between environments using the **conda activate** command
 * an environment may be **exported**, which simply means asking Conda to list all the packages in that environment, as well as the channels
   where they are available

Some key conda commands:

 * Make a new environment
 * List the environments currently known to conda
 * Switch to using an environment
 * Install a single package into the current environment
 * Print a list of all packages in the current environment in YAML format

Exercises -
Find out what version of fastx_toolkit is installed in the current conda environment.
Create a new conda env named ?? and install ??cutadapt?? into it.

### Using conda with Snakemake

Up to now, the `shell` commands you have been running with Snakemake are whatever is installed on the local machine. The commands may
have been installed with conda, or with the system package manager, or manually compiled. If you move your Snakefile to another machine,
or share it with a colleague, you will need to make sure all the dependencies are installed. If the versions of the packages are not the
same on two systems, you may discover that the same workflow produces different results.

If you are familiar with Conda, you may think to install the dependencies for your workflow into a conda environment, then `conda env export`
that into a YAML file, which you can use to quickly set up the same environment on another machine.

Snakemake takes this one step further with the following features:

1) Snakemake allows you to specify different environments to use for different rules. This is really useful if different rules need incompatible packages.
2) Snakemake can set up the environment for you, based on the YAML file

Export the env you built above
Then mention simplified YAML file

~~~
Simplified file here
~~~

Now let's make a very simple workflow that uses cutadapt on just one sequence. Yeah. That should do.
