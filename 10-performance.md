---
title: Optimising workflow performance
teaching: 20
exercises: 20
---

::::::::::::::::::::::::::::::::::::::: objectives

- Understand CPU, RAM and I/O bottlenecks
- Understand the *threads* declaration
- Use standard Linux tools to look at resource usage

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What compute resources are available on my system?
- How do I define jobs with more than one thread?
- How do I measure the compute resources being used by a workflow?
- How do I run my workflow steps in parallel?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the final Snakefile from episodes 1 to 8](files/ep08.Snakefile) you
may use to start this episode.*

## Processes, threads and processors

Some definitions:

- **Process** -     A running program (in our case, each Snakemake job can be considered
  one process)
- **Threads** -     Each process has one or more threads which run in parallel
- **Processor** -   Your computer has multiple *CPU cores* or processors, each of which can run one
  thread at a time

These definitions are a little simplified, but fine for our needs. The operating system kernel
shares out threads among processors:

- Having *fewer threads* than *processors* means you are not fully using all your CPU cores
- Having *more threads* than *processors* means threads have to "timeslice" on a core which is
  generally suboptimal

If you tell Snakemake how many threads each rule will use, and how many cores you have available,
it will start jobs in parallel to use all your cores. In the diagram below, five jobs are ready to
run and there are four system cores.

![][fig-threads]

## Listing the resources your Linux machine

Find out how many CPU cores you have on your machine with the `lscpu` command.

```bash
$ lscpu
```

Likewise find out the amount of RAM available:

```bash
$ free -h
```

And finally disk space, on the current partition:

```bash
$ df -h .
```

(or `df -h` without the `.` to show all partitions)

## Parallel jobs in Snakemake

You may want to see the relevant part of [the Snakemake documentation
](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads).

We'll force all the trimming and  kallisto steps to re-run by using the -F flag to Snakemake and
time the whole run using the standard `/usr/bin/time -v` command. You have to type the command like
this because `time` is a built-in command in BASH which takes precedence, so eg:

```bash
$ /usr/bin/time -v snakemake -j1 -F -- kallisto.{ref,temp33,etoh60}_{1,2,3}
```

:::::::::::::::::::::::::::::::::::::::  challenge

## Measuring how concurrency affects execution time

What is the *wallclock time* reported by the above command? We'll work out the average for the
whole class, or if you are working through the material on your own repeat the measurement three
times to get your own average.

Now change the Snakemake concurrency option to  `-j2` and then `-j4`.

- How does the total execution time change?
- What factors do you think limit the power of this setting to reduce the execution time?

:::::::::::::::  solution

## Solution

The time will vary depending on the system configuration but somewhere around 30 seconds is
expected, and this should reduce to around 25 secs with `-j2` but higher `-j` will produce
diminishing returns.

Things that may limit the effectiveness of parallel execution include:

- The number of processors in the machine
- The number of jobs in the DAG which are independent and can therefore be run in parallel
- The existence of single long-running jobs like *kallisto\_index*
- The amount of RAM in the machine
- The speed at which data can be read from and written to disk

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

There are **a few gotchas** to bear in mind when using parallel execution:

1. Parallel jobs will use more RAM. If you run out then either your OS will swap data to disk,
  or a process will crash
2. Parallel jobs may trip over each other if they try to write to the same filename at the same
  time (this can happen with temporary files, and in fact is a problem with our current
  `fastqc` rule definition)
3. The on-screen output from parallel jobs will be jumbled, so save any output to log files instead

## Multi-thread rules in Snakemake

In the diagram at the top, we showed jobs with 2 and 8 threads. These are defined by adding a
`threads:` part to the rule definition. We could do this for the *kallisto\_quant* rule:

```source
rule kallisto_quant:
    output:
        outdir = directory("kallisto.{sample}"),
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    threads: 4
    shell:
        "kallisto quant -t {threads} -i {input.index} -o {output.outdir} {input.fq1} {input.fq2}"
```

You should explicitly use `threads: 4` rather than `params: threads = "4"` because Snakemake
considers the number of threads when scheduling jobs. Also, if the number of threads requested for
a rule is less than the number of available processors then Snakemake will use the lower number.

We also added `-t {threads}` to the shell command. This only works for programs which allow you to
specify the number of threads as a command-line option, but this applies to a lot of different
bioinformatics tools.

:::::::::::::::::::::::::::::::::::::::  challenge

## Getting other programs to run with multiple threads

Find out how to set the number of threads for our *salmon\_quant* and *fastqc* jobs. Which of the
options below would need to be added to the shell command in each case?

1. `-t {threads}`
2. `-p {threads}`
3. `-num_threads {threads}`
4. multi-threaded mode is not supported

*Hint: use `salmon quant --help-alignment` and `fastqc --help`, or search the
online documentation.*

Make the corresponding changes to the Snakefile.

:::::::::::::::  solution

## Solution

For *salmon\_quant*, `-p {threads}` or equivalently `--threads {threads}` will work.

For *fastqc*, it may look like `-t {threads}` is good but this only sets "the number of files
which can be processed simultaneously", and the rule we have only processes a single file per
job. So in fact the answer is that, for our purposes, multi-threading is unsupported.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Fine-grained profiling

Rather than timing the entire workflow, we can ask Snakemake to benchmark an individual rule.

For example, to benchmark the `kallisto_quant` step we could add this to the rule definition:

```source
rule kallisto_quant:
    benchmark:
        "benchmarks/kallisto_quant.{sample}.txt"
    ...
```

The dataset here is so small that the numbers are tiny, but for real data this can be very
useful as it shows time, memory usage and IO load for all jobs.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Running jobs on a cluster

Learning about clusters is beyond the scope of this course, but for modern bioinformatics they are
an essential tool because many analysis jobs would take too long on a single computer. Learning to
run jobs on clusters normally means writing batch scripts and re-organising your code to be
cluster-aware. But if your workflow is written in Snakemake, it will run on a cluster will little
to no modification. Snakemake turns the individual jobs into cluster jobs, then submits and
monitors them for you.

To run Snakemake in this way, someone will need to determine the right parameters for your
particular cluster and save them as a profile. Once this is working, you can share the profile
with other users on the cluster, so discuss this with your cluster sysadmin.

Instructions for configuring the SLURM executor plugin can be found in the [Snakemake plugin
ctalog](https://snakemake.github.io/snakemake-plugin-catalog/), along with the *drmaa*,
*cluster-generic* and *cluster-sync* plugins which can support *PBS*, *SGE* and other cluster
schedulers.

[We have some specific suggestions for Eddie, the University of Edinburgh cluster](
files/snakemake_on_eddie.pdf)

![][fig-cluster]

::::::::::::::::::::::::::::::: instructor

## Running on cluster and cloud

Running workflows on HPC or Cloud systems could be a whole course in itself. The topic is too
important not to be mentioned here, but also complex to teach because you need a cluster to work
on.

The original author of this material would demonstrate running the workflow on the [Cirrus](
https://cirrus.ac.uk) system in Edinburgh. If you are teaching this lesson and have institutional
HPC then ideally you should liaise with the administrators of the system to make a suitable
installation of a recent Snakemake version and a profile to run jobs on the cluster job scheduler.
In practise this may be easier said than done!

If you are able to demonstrate Snakemake running on cloud as part of one of these courses then
we'd much appreciate any feedback on how you did this and how it went.

::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Cluster demo

A this point in the course there may be a cluster demo...

::::::::::::::::::::::::::::::::::::::::::::::::::

[comment]: # (
Photo credit: Cskiran
Sourced from Wikimedia Commons
CC-BY-SA-4.0
)

*For reference, [this is a Snakefile](files/ep10.Snakefile) incorporating the changes made in
this episode.*



[fig-threads]: fig/snake_threads.svg {alt='Representation of a computer with four microchip icons
indicating four available cores. To the right are five small green boxes representing Snakemake
jobs and labelled as wanting 1, 1, 1, 2 and 8 threads respecively.'}
[fig-cluster]: fig/Multiple_Server_.jpg {alt='A photo of some high performance computer hardware
racked in five cabinets in a server room. Each cabinet is about 2.2 metres high and 0.8m wide. The
doors of the cabinets are open to show the systems inside. Orange and yellow cabling is
prominent, connecting ports within the second and third racks.'}


:::::::::::::::::::::::::::::::::::::::: keypoints

- To make your workflow run as fast as possible, try to match the number of threads to the number of cores you have
- You also need to consider RAM, disk, and network bottlenecks
- Profile your jobs to see what is taking most resources
- Snakemake is great for running workflows on compute clusters

::::::::::::::::::::::::::::::::::::::::::::::::::


