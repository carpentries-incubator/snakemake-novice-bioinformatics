---
title: "Robust quoting in Snakefiles"
teaching: 0
exercises: 0
questions:
- "How are shell commands processed before being run?"
- "How do I avoid quoting problems in filenames and commands?"
objectives:
- "Review quoting rules in the shell"
- "Understand how shell command strings are processed"
- "Understand how to make Snakemake commands robust"
keypoints:
- "Having a handle on string quoting is boring but important"
- "Understand the three processing steps each *shell* command goes through before it is actually run"
- "Make use of triple-quotes"
- "Watch out for commands that have {curly brackets} and double them up"
- "Use the built-in `:q` feature to protect arguments from Bash interpretation"
---

## A review of quoting rules in the Bash shell

Consider the following simple Bash shell command:

~~~
$ echo Why is a "mouse" when     it spins?
Why is a mouse when it spins?
~~~

The message is printed back, but not before the shell has interpreted the text as per the standard command-line
parsing rules

 * The quotes around '"mouse"' have been removed
 * The extra spaces between 'when' and 'if' have been ignored
 * More subtly, the last word will be interpreted as a glob pattern...

~~~
$ touch spinsX spinsY
$ echo Why is a "mouse" when     it spins?
Why is a mouse when it spinsX spinsY
~~~

*Note: if you have certain shell settings you may have seen a warning about the unmatched glob pattern.*

In shell commands, some characters are "safe", in that Bash does not try to interpret them at all:

* Letters `A-Z` and `a-z`
* Digits `0-9`
* Period, underscore, forward slash, and hyphen `._-`

Most other characters have some sort of special meaning and must be enclosed in quotes if you want to pass them through
verbatim to the program you are running. This is essential with things like `awk` commands:

~~~
# Print mean read length in FASTQ file - from https://github.com/stephenturner/oneliners
$ awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' reads/ref_1_1.fq
101
~~~

The 'single quotes' ensure that the expression is passed directly to `awk`. If double quotes had been used, Bash would
replace "$0" with the name of the current script and break the awk logic, giving a meaningless result.

~~~
# With double quotes the "$0" is substituted by Bash, rather than being passed on to awk.
$ awk "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}" reads/ref_1_1.fq
0
# In an interactive shell $0 is a variable containing "bash"
$ echo "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}"
NR%4==2{sum+=length(bash)}END{print sum/(NR/4)}
~~~

So in Bash, putting a string in single quotes is normally enough to preserve the contents, but if you need to add literal
single quotes to the awk command for any reason you have to do some awkward (pun intended) construction.

~~~
# Yuk! But it does work...
$ awk '{print "\"double\" '"'"'single'"'"'"}' <<<''
"double" 'single'
~~~

## Quoting rules in Snakemake

Anyone who has done any amount of programming or shell scripting will have come across quoting issues in code. Bugs related
to quoting can be very troublesome. In Snakemake these are particularly complex because the "shell" part of each rule
undergoes three rounds of interpretation before anything is actually run:

 1. The string is parsed according to Python quoting rules
 1. Placeholders in curly brackets (eg. {input} {output}) are then replaced
 1. The resulting string goes to the Bash shell and is subject to all Bash parsing rules

We'll now look at some best practises for making your Snakefiles robust, and some simple rules to avoid most
mis-quoting complications.

> ## Exercise - adding a lenreads rule
>
> Say we add a new rule named *lenreads*, looking very much like the existing *countreads* rule and using the awk
> expression we saw earlier.
>
> ~~~
> rule lenreads:
>   output: "{indir}.{asample}.fq.len"
>   input:  "{indir}/{asample}.fq"
>   shell:
>       "awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' {input} > {output}"
> ~~~
>
> Will this work as shown? If not, why not? Try it and see.
>
> > ## Answer
> >
> > It won't work. Snakemake assumes that all parts of the string in {curlies} are placeholders.  The error will say
> > something like `NameError: The name 'sum+=length($0)' is unknown in this context`. To resolve this, we
> > have to double up all the curly braces:
> >
> > ~~~
> >   shell:
> >     "awk 'NR%4==2{{"{{"}}sum+=length($0)}}END{{"{{"}}print sum/(NR/4)}}' {input} > {output}"
> > ~~~
> {: .solution}
{: .challenge}

## Best practise for quoting

For complex commands, use the *triple-quoted r-strings* we saw earlier.

~~~
r"""Strings like this"""
~~~

The syntax allows embedded newlines, literal \n \t, and both types of '"quotes"'. In other words, the interpretation
as a Python string does as little as possible, leaving most interpretation to the shell.

The triple-quoting does not protect {curlies}, so if you are needing to use awk commands like the one above, rather
the adding extra braces into the command you can put it in a variable.

~~~
LEN_READS_CMD = "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}"

rule lenreads:
  shell:
      "awk '{LEN_READS_CMD}' {input} > {output}"
~~~

Or even better:

~~~
rule lenreads:
  shell:
      "awk {LEN_READS_CMD:q} {input} > {output}"
~~~

Using `{LEN_READS_CMD:q}` instead of `'{LEN_READS_CMD}'` is asking Snakemake to quote the awk command for you. In
this case, Snakemake will just put it into single quites, but if your variable contains single quotes or embedded
newlines or tabs or any other oddities then Snakemake will quote it robustly for you.

The `:q` syntax works on any placeholder and you can safely add it to all the placeholders, so we could well say:

~~~
rule lenreads:
  shell:
      "awk {LEN_READS_CMD:q} {input:q} > {output:q}"
~~~

Now the *lenreads* rule would be able to work on an input file that contains spaces or other unusual characters,
and if the *input* is a list of files this will still work just fine, whereas `'{input}'` will fail as it just
combines all the filenames into one big string.

In general, choose file names that only contain shell-safe characters and no spaces, but if you can't do that
then just ensure all your placeholders have `:q` and you should be fine.

> ## Exercise
>
> Use the following command to rename the *temp33* and *etoh60* samples so that the filenames contain spaces:
>
> ~~~
> $ rename -v -s 'etoh' 'etoh ' -s 'temp' 'temp ' reads/*.fq
> ~~~
>
> Fix the workflow so you can still run the MultiQC report over all the samples (run Snakemake with `-F` to check
> that all the steps work).
>
> > ## Solution
> >
> > This just involves adding `:q` to a whole bunch of placeholders. Unless you are very diligent it will probably
> > take a few tries to catch every one of them and make the workflow functional again.
> >
> {: .solution}
{: .challenge}

{% include links.md %}

