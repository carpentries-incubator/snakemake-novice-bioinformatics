---
title: "Robust quoting in Snakefiles"
teaching: 20
exercises: 20
questions:
- "How are shell commands processed before being run?"
- "How do I avoid quoting problems in filenames and commands?"
objectives:
- "Review quoting rules in the shell"
- "Understand how shell command strings are processed"
- "Understand how to make Snakemake commands robust"
keypoints:
- "Having a grasp of string quoting rules is boring but important"
- "Understand the three processing steps each *shell* command goes through before it is actually run"
- "Make use of triple-quotes"
- "Watch out for commands that have {curly brackets}, and double them up"
- "Use the built-in `:q` feature to protect arguments from *bash* interpretation"
---
*For reference, [this is the final Snakefile from episodes 1 to 7](../code/ep07.Snakefile) you may use to
start this episode.*

## A review of quoting rules in the *bash* shell

Consider the following simple *bash* shell command:

~~~
$ echo Why is a "mouse" when     it spins?
Why is a mouse when it spins?
~~~
{: .language-bash}

The message is printed back, but not before the shell has interpreted the text as per various command-line
parsing rules

 * The quotes around `"mouse"` have been removed
 * The extra spaces between `when` and `it` have been ignored
 * More subtly, `spins?` will be interpreted as a glob pattern if there is any possible match...

~~~
$ touch spinsX spinsY
$ echo Why is a "mouse" when     it spins?
Why is a mouse when it spinsX spinsY
~~~
{: .language-bash}

*Note: if you have certain shell settings you may have seen a warning about the unmatched glob pattern.*

In shell commands, some characters are "safe", in that *bash* does not try to interpret them at all:

* Letters `A-Z` and `a-z`
* Digits `0-9`
* Period, underscore, forward slash, and hyphen `._/-`

Most other characters have some sort of special meaning and must be enclosed in quotes if you want to pass them through
verbatim to the program you are running. This is essential with things like *awk* commands:

~~~
# Print mean read length in FASTQ file - from https://github.com/stephenturner/oneliners
$ awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' reads/ref_1_1.fq
101
~~~
{: .language-bash}

The 'single quotes' ensure that the expression is passed directly to *awk*. If "double quotes" had been used, *bash* would
replace `$0` with the name of the current script and break the awk logic, giving a meaningless result.

~~~
# With double quotes the "$0" is substituted by *bash*, rather than being passed on to awk.
$ awk "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}" reads/ref_1_1.fq
0
# In an interactive shell $0 is a variable containing the value "bash"
$ echo "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}"
NR%4==2{sum+=length(bash)}END{print sum/(NR/4)}
~~~
{: .language-bash}

So in *bash*, putting a string in single quotes is normally enough to preserve the contents, but if you need to add literal
single quotes to the awk command for any reason you have to do some awkward (pun intended) construction.

~~~
# Yuk! But it does work...
$ awk '{print "\"double\" '"'"'single'"'"'"}' <<<''
"double" 'single'
~~~
{: .language-bash}

## Quoting rules in Snakemake

Anyone who has done any amount of programming or shell scripting will have come across quoting issues in code. Bugs related
to quoting can be very troublesome. In Snakemake these are particularly complex because the "shell" part of each rule
undergoes three rounds of interpretation before anything is actually run:

 1. The string is parsed according to Python quoting rules
 1. Placeholders in curly brackets (eg. `{input}` `{output}`) are then replaced
 1. The resulting string goes to the *bash* shell and is subject to all *bash* parsing rules

We'll now look at some best practises for making your Snakefiles robust, and some simple rules to avoid most
mis-quoting complications.

> ## Exercise - adding a *lenreads* rule
>
> Say we add a new rule named *lenreads*, looking very much like the existing *countreads* rule and using the *awk*
> expression we saw earlier.
>
> ~~~
> rule lenreads:
>     output: "{indir}.{myfile}.fq.len"
>     input:  "{indir}/{myfile}.fq"
>     shell:
>         "awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' {input} > {output}"
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
> >       "awk 'NR%4==2{{"{{"}}sum+=length($0)}}END{{"{{"}}print sum/(NR/4)}}' {input} > {output}"
> > ~~~
> {: .solution}
{: .challenge}

## Best practise for writing robust workflows

For complex commands in Snakemake rules, use the *triple-quoted r-strings* we saw earlier.

~~~
r"""Strings like this"""
~~~

The syntax allows embedded newlines, literal \n \t, and both types of '"quotes"'. In other words, the interpretation
as a Python string does as little as possible, leaving most interpretation to the *bash* shell.

The triple-quoting does not protect {curlies}, so if you are needing to use *awk* commands like the one above, rather
than adding extra braces into the command you could define it as a variable.

~~~
LEN_READS_CMD = r"""NR%4==2{sum+=length($0)}END{print sum/(NR/4)}"""

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
this case, Snakemake will just put it into single quotes, but if your variable contains single quotes or embedded
newlines or tabs or any other oddities then Snakemake will quote it robustly for you.

The `:q` syntax works on any placeholder and you can safely add it to all the placeholders, so we could well say:

~~~
rule lenreads:
    shell:
        "awk {LEN_READS_CMD:q} {input:q} > {output:q}"
~~~

Now the *lenreads* rule would be able to work on an input file that contains spaces or other unusual characters.
Also, if the *input* is a list of files, this will still work just fine, whereas `'{input}'` will fail as it just
combines all the filenames into one big string within the quotes.

In general, choose file names that only contain shell-safe characters and no spaces, but if you can't do that
then just ensure all your placeholders have `:q` and you should be fine.

> ## Exercise
>
> Use the following command to rename the *temp33* and *etoh60* samples so that the filenames contain spaces:
>
> ~~~
> $ rename -v -s 'etoh' 'etoh ' -s 'temp' 'temp ' reads/*.fq
> ~~~
> {: .language-bash}
>
> Fix the workflow so you can still run the MultiQC report over all the samples (run Snakemake with `-F` to check
> that all the steps really work).
>
> > ## Solution
> >
> > This just involves adding `:q` to a whole bunch of placeholders. Unless you are very diligent it will probably
> > take a few tries to catch every one of them and make the workflow functional again.
> >
> {: .solution}
{: .challenge}

*For reference, [this is a Snakefile](../code/ep13.Snakefile) incorporating the changes made in this episode.*

{% include links.md %}

