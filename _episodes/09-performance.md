---
title: "Optimising workflow performance"
teaching: 20
exercises: 10
questions:
- "What compute resources are available on my system?"
- "How do I measure the compute resources being used by a workflow?"
- "How do I run my workflow steps in parallel?"
- "..."
objectives:
- "Understand CPU, RAM and I/O bottlenecks"
- "Use standard Linux tools to look at resource usage"
keypoints:
- "Add key points"
---

## Processes, threads and processors

Some definitions:

* **Process** - 	A running program (in our case, each Snakemake job can be considered one process)
* **Threads** - 	Each process has one or more threads which run in parallel
* **Processor** -	Your computer has multiple *CPU cores* or processors, each of which can run one thread at a time

These definitions are a little simplified, but fine for our needs. The operating system kernel shares out threads among processors:

* Having *fewer threads* than *processors* means you are not fully using all your CPU cores
* Having *more threads* than *processors* means threads have to share a core which is generally suboptimal

If you tell Snakemake how many threads each rule will use, and how many cores you have available, it will start jobs
in parallel to use all your cores.

![Allocating cores to jobs in Snakemake][fig-threads]

## Profiling your Linux machine

Find out how many CPU cores you have on your machine with the `lscpu` command.

~~~
$ lscpu
~~~

Likewise find out the amount of RAM available:

~~~
$ free -h
~~~

And finally disk space, on the current partition:

~~~
$ df -h .
~~~

(or `df -h` to show all partitions)

## Parallel jobs in Snakemake

You may want to see the relevant part of
[the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads).

We'll force all the trimming and  kallisto steps to re-run by using the -F flag to Snakemake and time
the whole run using the standard `/usr/bin/time -v` command. You have to type the command like this because
`time` is a built-in command in BASH which takes precedence, so eg:

~~~
$ /usr/bin/time -v snakemake -j1 -F -- kallisto.{ref,temp33,etoh60}_{1,2,3}
~~~


> ## Exercise
>
> What is the 'wallclock time' reported by the above command? We'll work out the average for the whole class, or
> if you are working through the material on your own repeat the measurement three times to get your own average.
>
> Now change the Snakemake concurrency option to  `-j 2` and then `-j 4`. How does the total execution time change?
> What do you think limits the power of this setting to reduce the execution time?
>
> > ## Solution
> >
> > The time will vary depending on the system configuration but somewhere around 30 seconds is expected, and this
> > should reduce to around 25 secs with `-j 2` but higher `-j` will produce diminishing returns.
> >
> > Things that may limit the effectiveness of parallel execution include:
> >
> > * The number of processors in the machine
> > * The number of jobs in the DAG which are independent and can therefore be run in parallel
> > * The existence of single long-running jobs like *kallisto_index*
> > * The amount of RAM in the machine
> > * The speed at which data can be read from and written to disk
> >
> {: .solution}
{: .challenge}


> ## Fine-grained profiling
>
> Rather than timing the entire workflow, we can ask Snakemake to benchmark an individual rule.
>
> For example, to benchmark the `kallisto_quant` step we could add this to the rule definition:
>
> ~~~
> benchmark:
> 	"benchmarks/kallisto_quant.{sample}.txt"
> ~~~
>
> The dataset here is so small that the numbers are tiny, but for real data this can be very useful as it shows time, memory
> usage and IO load for all jobs.
>
>
{: .callout}

There are **a few gotchas** to bear in mind when using parallel execution:

1. Parallel jobs will use more RAM. If you run out then either your OS will swap data to disk, or a process will crash
1. Parallel jobs may trip over each other if they try to write to the same file at the same time (can happen with temporary files)
1. The on-screen output from paralle jobs will be jumbled, so save output to log files instead

## Running jobs on a cluster

Learning about clusters is beyond the scope of this course, but for modern bioinformatics they are an essential tool because
many analysis jobs would take too long on a single computer. Learning to run jobs on clusters normally means writing batch
scripts and re-organising your code to be cluster-aware. But if your workflow is written in Snakemake, it will run on a cluster
will little to no modification. Snakemake turns the jobs into cluster jobs, then submits and monitors them for you.

![Some high performance compute][fig-cluster]


> ## Cluster demo
>
> A this point in the course there may be a cluster demo...
>
{: .callout}

[fig-threads]: ../fig/snake_threads.svg
[fig-cluster]: ../fig/Multiple_Server_.jpg

{% include links.md %}
