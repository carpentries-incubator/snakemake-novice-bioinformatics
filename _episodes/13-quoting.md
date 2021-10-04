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
- "Understand issues regarding quoting special characters"
keypoints:
- "Add key points"
---

## A review of quoting rules in the Bash shell

Consider the following simple Bash shell command:

~~~
$ echo How "cool" are    you?
How cool are you?
~~~

The message is printed back, but not before the shell has interpreted the text as per the standard command-line
parsing rules

 * The quotes around '"cool"' have been removed
 * The extra spaces between 'are' and 'you' have been ignored
 * More subtly, the last word will be interpreted as a glob pattern...

~~~
$ touch youX youY
$ echo How "cool" are    you?
How cool are youX youY
~~~

*Note: if you have certain shell settings you may see a warning about the unmatched glob pattern.*

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
# With double quotes the "$0" is replaced by Bash before awk can see it.
$ awk "NR%4==2{sum+=length($0)}END{print sum/(NR/4)}" reads/ref_1_1.fq
0
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
 1. The resulting string goes to the shell and is subject to all Bash parsing rules

We'll now look at some best practises for making your Snakefiles robust, and some simple rules to avoid most
mis-quoting complications.

###### Now have some exercises that illustrate the three points and introduce ways to deal with them:

> ## Exercise - adding a lenreads rule
>
> Say we add a new rule named 'lenreads', looking very much like the existing 'countreads' rule and using the awk
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
> Will this work as shown? If not, why not?
>
> > ## Answer
> >
> > It won't work. Snakemake assumes that all parts of the string in {curlies} are placeholders. In this case, we
> > have to double them up:
> >
> > ~~~
> >   shell:
> >     "awk 'NR%4==2{{sum+=length($0)}}END{{print sum/(NR/4)}}' {input} > {output}"
> > ~~~

## Best practise for quoting

For complex commands, use the triple-quote strings we saw earlier.

~~~
r"""Strings like this"""
~~~

The syntax allows embedded newlines, literal \n \t, and both types of '"quotes"'. In other words, the interpretation
as a Python string does as little as possible, leaving most interpretation to the shell.

Choose file names that only contain shell-safe characters and no spaces. This is good general proactise in any case.


Double-escaping {{curlies}}. There is no good way to avoid this. Or is there? We can add them outside of the rule, maybe?
The AWK example is good here.

~~~
LEN_READS_CMD = ['awk', 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}']

rule lenreads:
  output: "{indir}.{asample}.fq.len"
  input:  "{indir}/{asample}.fq"
  shell:
      "{LEN_READS_CMD:q} {input} > {output}"

~~~

Shell parsing rules. Normally you want this because you want your pipes, redirects, shell variables, etc. to work.
If you choose filenames using only [0-9A-Z-a-z.\_] then you will always be "safe", but sometimes you have to deal with
awkwardly named files or parameters. If you want them unmolested. Use {foo:q}. Not '{foo}'!
Have an exercise that demonstrates why this is good.

Getting increasingly unsure about this. Nead more hummus to work out how to make this a good chapter.

How will I best demonstrate the power of {:q}?

letters = ['"', "'", "$ ", "\n\tx", "{:}"]
shell:
    "echo {letters:q}"

Hmmm. A little obscure. Well have a think.


{% include links.md %}

