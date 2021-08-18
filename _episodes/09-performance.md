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

Soe definitions:

* **Process** - 	A running program (in our case, a Snakemake job)
* **Threads** - 	Each process has one or more threads
* **Processor** -	Your computer has multiple CPU cores, which may run threads or processes in parallel

Linux shares out threads among processors:

* Having fewer threads than processors means you are not fully using all your cores
* Having more threads than processors is generally suboptimal

If you tell Snakemake how many threads a rule will use, and how many cores you have available, it will start jobs
in parallel to use all your cores.

![Allocating cores to jobs in Snakemake][fig-threads]

## Profiling your Linux machine

Find out how many CPU cores you have on your VM with the `lscpu` command.

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

## Performance and Resource Usage

You may want to see the relevant part of
[the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads).


Optimising CPU Usage

We'll force all the trimming and  kallisto steps to re-run by using the -F flag to Snakemake and time
the whole run using the standard `/usr/bin/time -v` command. You have to type the command like this because
`time` is a built-in command in BASH which takes precedence, so eg:

~~~
$ /usr/bin/time -v snakemake -j1 -F -- kallisto.{ref,temp33,etoh60}_{1,2,3}
~~~

> ## Exercise
> Do this two more times. What is the average 'wallclock time' over the three?
> Now try again, using the Snakemake option "-j 2". How does the execution time change? What about "-j 3" or "-j 4".
> What limits the power of this setting to reduce the execution time?
>
{: .challenge}


Rather than timing the entire workflow, we can ask Snakemake to benchmark just one rule.

To benchmark the `kallisto_quant` step we add this to the rule definition:

~~~
benchmark:
	"benchmarks/kallisto_quant.{sample}.txt"
~~~

How long is Kallisto taking to run? (For simplicity, just use "-j 1" for all these tests.)

Now increase the number of threads used by BWA to 2 by setting the "-t 2" option in the "bwa aln" command. What impact does this have?

FIXME - Kallisto runs too quickly :-(. I need to find a better example here.

What limits our ability to speed up workflows in this way?

## Optimising memory usage

Monitoring memory use is actually quite tricky. The /usr/bin/time -v command gives you a "Maximum resident set size" which is the amount
of RAM used, but when you timed the whole pipeline only Snakemake itself (ie. the Python runtime) was actually monitored. Most of the
resources were being used by child processes like `kallisto` and were never counted. Snakemake benchmarking can be used to profile the
individual tasks, but to find the overall memory usage of your workflow a simpler option is to monitor the total free memory on the system
and then see how much is used up as the pipeline runs. Various memory monitoring tools are available for Linux, but here is a simple method
that works with just the basic tools on the Linux system and a spreadsheet to visualise the result.

Exercise G.4:

In a separate terminal window, run this monitoring command and leave it running:

~~~
$ ( while true ; do grep 'MemAvailable' /proc/meminfo ; \
  sleep 0.5 ; done ) | tee ~/memavail.log
~~~

Now while the monitor is running, re-run Snakemake with the "-j 1" option and when that finishes run it with the "-j 3" option.

Now stop the memory monitoring (with Ctrl+C).

Open the LibreOffice Calc spreadsheet from the Applications menu and open the memavail.log file. You'll need to set the delimiter to
"space" so that all the numbers end up in one column. Select just this column and make a chart (Insert > Chart...). Make it a line graph.

How does the shape of the graph reflect the processes being run by Snakemake? How does the memory usage compare between the two runs?
Why do we see this shape?

### FIXME - this doesn't work too well with the example we have. Work it out.

[fig-threads]: ../fig/snake_threads.svg

{% include links.md %}
