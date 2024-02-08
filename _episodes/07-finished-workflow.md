---
title: "Finishing the basic workflow"
teaching: 40
exercises: 40
questions:
- "What does the full sample workflow look like?"
- "How can we report some initial results from this analysis?"
objectives:
- "Finish the sample workflow to produce a final MultiQC report"
- "See more ways to handle awkward software by adding extra setup shell commands"
keypoints:
- "Once a workflow is complete and working, there will still be room for refinement"
- "This completes the introduction to the fundamentals of Snakemake"
---
*For reference, [this is the Snakefile](../code/ep06.Snakefile) you should have to start
the episode.*

We've seen how to link rules in a pipeline and how to merge all the results at the final step. This
is the basic pattern for many analysis workflows. For simplicity, in episode 5, we just used the
`cat` command to combine all the `.count` files but now we'll use *MultiQC* to combine the results
from Kallisto and FastQC into a single report for all samples. We'll also add in an alternative
quantification tool called *Salmon* and this will complete the pipeline.

## The full workflow we are constructing

* **fastq\_quality\_trimmer** is part of the FastX toolkit and removes low-quality basecalls from
  the raw reads.
  *We first used it in [episode 2]({{ page.root }}{% link _episodes/02-placeholders.md %}).*
* **FastQC** calculates a variety of metrics on a FASTQ file and produces an HTML report and a
  ZIP file.
  *We introduced this in [episode 6]({{ page.root }}{% link _episodes/06-awkward_programs.md %}).*
* **Kallisto** aligns the reads to a reference transcriptome and produces a table of
  transcript abundance.
  *We first used it in [episode 3]({{ page.root }}{% link _episodes/03-chaining_rules.md %}).*
* **Salmon** is a alternative to Kallisto, using a different alignment algorithm.
  *We've not used it yet.*
* **MultiQC** combines the reports from various tools, including FastQC, Kallisto, and Salmon,
  into a single HTML report over all samples. This is by no means a full RNA-Seq analysis report
  but today it completes our Snakemake pipeline.

![Summary of our full QC workflow with icons representing the steps listed above.
The input data is also summarized, with 18 paired FASTQ files under yeast/reads, for
the three repeats of all three conditions, as well as the transcriptome in gzipped
FASTA format.
][fig-workflow]

At this point we have everything we need, in terms of Snakemake knowledge, to add the two remaining
tools and complete the Snakefile. As with FastQC, some quirks in the way the new tools work will
need to be accounted for.

## Adding Salmon as an alternative to Kallisto

At the end of the previous episode we modified the *kallisto* rule by declaring that the output of
the rule was a directory, rather than explicitly listing the three files that Kallisto writes in
that directory. You should ensure that you are working with this new version of the *kallisto* rule
as it will be a template for adding a *salmon* rule.

> ## Exercise - adding Salmon as an alternative to Kallisto
>
> An alternative application for transcript quantification is Salmon. The procedure is virtually
> identical, to Kallisto, having an indexing step and a quantification step.
> Note that in real usage you are advised to prepare and add decoy sequences to the transcriptome
> index, but for the purposes of this tutorial we'll just keep things as simple as possible.
>
> Based upon the following command templates:
>
> ~~~
> $ salmon index -t <transcriptome as fasta> -i <index name> -k 31
> $ salmon quant -i <index name> -l A -1 <fastq1> -2 <fastq2> --validateMappings -o <output path>
> ~~~
> {: .language-bash}
>
> Add a pair of rules to index and align the reads with Salmon. Note that:
>
> 1. Unlike Kallisto, the index produced by Salmon is a directory of files, not a single file - so
>    both of these new rules will produce a directory as output.
> 1. As noted in the last episode, you only need the `directory()` marker on the outputs of rules,
>    not the inputs.
>
> > ## Solution
> >
> > ~~~
> > rule salmon_quant:
> >     output: directory("salmon.{sample}")
> >     input:
> >         index = "Saccharomyces_cerevisiae.R64-1-1.salmon_index",
> >         fq1   = "trimmed/{sample}_1.fq",
> >         fq2   = "trimmed/{sample}_2.fq",
> >     shell:
> >         "salmon quant -i {input.index} -l A -1 {input.fq1} -2 {input.fq2} --validateMappings -o {output}"
> >
> > rule salmon_index:
> >     output:
> >         idx = directory("{strain}.salmon_index")
> >     input:
> >         fasta = "transcriptome/{strain}.cdna.all.fa.gz"
> >     shell:
> >         "salmon index -t {input.fasta} -i {output.idx} -k 31"
> > ~~~
> > {: .source}
> >
> > If you copied the *kallisto_index* rule and logged the output of *salmon_index* to a file this
> > is fine. Just make sure when copying and pasting that you change all the parts that need
> > to change!
> >
> {: .solution}
{: .challenge}


## Combining the outputs with MultiQC

MultiQC scans for analysis report files in a given directory and all subdirectories, then makes a
report of everything it finds. It knows about FastQC, Salmon and Kallisto outputs so we should be
able to compile a report on all these. To try this out and scan the current directory, simply run:

~~~
$ multiqc . -o multiqc_out
~~~
{: .language-bash}

> ## Exercise - adding a MultiQC rule
>
> Earlier, in episode 5, we made a basic summary-type rule called *all_counts*. Now make a
> *multiqc* rule that gathers up all the FastQC, Salmon and Kallisto reports.
>
> Considerations:
>
> 1. Your rule is going to have several named inputs, and these inputs will be lists of files
>    generated with `expand()` functions.
> 1. Ensure that both *kallisto_quant* and *salmon_quant* are run on all 9 samples, that is all
>    three repeats of all three conditions.
> 1. Run FastQC on the untrimmed reads only, and note that MultiQC specifically uses the `.zip`
>    files for input, not the `.html`.
> 1. Since multiqc scans for input files, the input names don't have to be explicitly mentioned in
>    the `shell` part.
>
> > ## Solution
> >
> > ~~~
> > rule multiqc:
> >     output: directory('multiqc_out')
> >     input:
> >         salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
> >     shell:
> >         "multiqc . -o multiqc_out"
> > ~~~
> > {: .source}
> >
> > Since the rule has no wildcards, you can run it by giving either the rule name or the output
> > directory name as a target.
> >
> > ~~~
> > $ snakemake -j1 -p multiqc
> >
> > ...or equivalently...
> >
> > $ snakemake -j1 -p multiqc_out
> > ~~~
> > {: .language-bash}
> >
> {: .solution}
{: .challenge}

> ## Exercise - fixing Kallisto
>
> You may notice that MultiQC is not capturing any Kallisto output when making the reports. The
> reason for this is given in the [MultiQC manual here](https://multiqc.info/docs/#kallisto):
>
> > *Note - MultiQC parses the standard out from Kallisto, not any of its output files
> > (abundance.h5, abundance.tsv, and run_info.json). As such, you must capture the Kallisto
> > output to a file when running it for use with MultiQC.*
> {: .callout}
>
>
> Fix the Snakefile so that Kallisto terminal output is redirected to a file and can be collected
> by MultiQC.
>
>  * *Hint 1:* The manual above is not quite right - you need to capture both
>    **stdout and stderr**, so use `>&` rather than `>`, as we did with the indexing step.
>  * *Hint 2:* MultiQC does not mind what you call the file, so choose your own sensible name.
>
> > ## Solution
> >
> > ~~~
> > # Kallisto quantification of one sample, with log capture.
> > rule kallisto_quant:
> >     output: directory("kallisto.{sample}")
> >     input:
> >         index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
> >         fq1   = "trimmed/{sample}_1.fq",
> >         fq2   = "trimmed/{sample}_2.fq",
> >     shell:
> >         """mkdir {output}
> >            kallisto quant -i {input.index} -o {output} {input.fq1} {input.fq2} >& {output}/kallisto_quant.log
> >         """
> > ~~~
> > {: .source}
> >
> > There are several perfectly good ways of structuring this, so don't worry if your answer is
> > different.
> >
> > A gotcha with the above version is that the output directory needs to be created before
> > *kallisto quant* is run, much like with FastQC.
> > Remember that Snakemake deletes any existing outputs, including outputs that are directories,
> > before the job is run, and while Kallisto will create the directory for you this is too late
> > for the shell to make the log file and without the `mkdir` command it will report
> > "No such file or directory".
> >
> > Another option is to make the log file outside of the directory, or stick with declaring the
> > individual files as outputs, as when we first made the rule, in which case the directory will
> > be made by Snakemake. It's also possible to declare both the directory and a file within that
> > directory as separate outputs, which is unusual but may be a good approach here.
> >
> {: .solution}
{: .challenge}

> ## Exercise - making the MultiQC rule more robust
>
> Because MultiQC scans for suitable input rather than taking an explicit list of files, there is a
> risk that it picks up unwanted reports if you leave old files sitting in the directory. To make
> the rule fully explicit, one idea is to make a temporary directory and symlink all the files into
> it, and then tell MultiQC to look in there. Amend the *multiqc* rule so it does this.
>
>  * *Hint 1:* When making links, use `ln -snr -t <target_dir> <src>` to avoid link
>    relativity issues.
>  * *Hint 2:* If you feel that tweaking some other rules would make this easier, feel free to
>    do that.
>
> > ## Solution
> >
> > This solution will work with the version of *kallisto_quant* in the solution above.
> >
> > ~~~
> > rule multiqc:
> >     output:
> >         mqc_out = directory('multiqc_out'),
> >         mqc_in  = directory('multiqc_in'),
> >     input:
> >         salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
> >         fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
> >     shell:
> >         """mkdir {output.mqc_in}
> >            ln -snr -t {output.mqc_in} {input}
> >            multiqc {output.mqc_in} -o {output.mqc_out}
> >         """
> > ~~~
> > {: .source}
> {: .solution}
{: .challenge}

To see the actual MultiQC report, open the file *multiqc_out/multiqc_report.html* in a web browser.
You can do this directly from the command line, assuming you have a default browser configured.

On Linux environments:

~~~
$ xdg-open multiqc_out/multiqc_report.html
~~~
{: .language-bash}

For MacOS:

~~~
$ open multiqc_out/multiqc_report.html
~~~
{: .language-bash}


The report has a few issues, but we'll not get distracted by the details of how to configure
MultiQC to resolve them.

[fig-workflow]: ../fig/overview.svg

*For reference, [this is a Snakefile](../code/ep07.Snakefile) incorporating the changes made in
this episode. You may now proceed to any later episode in the lesson using this workflow as a
starting point.*

{% include links.md %}
