---
title: 'Reference'
---

## Glossary

- **Snakemake**

A free software application, written in [Python](https://swcarpentry.github.io/python-novice-inflammation/),
which knows how to read and execute Snakefiles. Once installed,
it is invoked from the shell via the *snakemake* command.

- **Snakefile**

A text file containing one or more rule definitions which may be invoked to make a workflow. The file may also
hold variables, functions and other settings.

- **Workflow**

A series of jobs that are to be executed in a suitable order to generate a final target output. Snakemake plans the
workflow based upon the Snakefile, the configuration, and the target, so a single Snakefile may define multiple
possible workflows.

- **DAG**

Directed Acyclic Graph. An internal representation of the jobs in a workflow and how they depend on one another. It
may be visualised as a diagram of boxes connected by arrows.

- **Rule**

A recipe telling Snakemake how to generate some particular type of file, normally by running a shell command on one or
more input files.

- **Job**

A concrete task in a workflow, also called a **Step** in the Snakemake logs. A job is defined by a rule plus specific
wildcard values, so each rule may yield many jobs.

- **Conda**

A package installation system derived from the [Anaconda distribution](https://www.anaconda.com/products/distribution),
which allows for management of self-contained software environments
and backed by a large repository of [community maintained packages](https://anaconda.org/search).

- **Cluster**

A computer system comprised of multiple processing nodes, connected to a shared storage system and supporting the running
of batch jobs via a cluster manager. When run with appropriate settings on a cluster, Snakemake can turn its own jobs
into cluster batch jobs and submit and monitor them automatically.

- **YAML**

A way to express structured data within a text file. Snakemake uses this format for configuration files, and Conda uses
it for exporting and importing environment definitions.

- **Escaping**

With regard to the text within a Snakefile, a way to indicate that certain characters should not be treated as having a
special meaning. A typical example would be `"\t"` versus `"\\t"`, where the first string represents a tab character,
but the second represents a backslash character followed by a letter 't', and the extra backslash is referred to as
an *escape character*.




