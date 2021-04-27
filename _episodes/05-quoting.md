---
title: "Quoting"
teaching: 0
exercises: 0
questions:
- "How do I make my Snakefiles robust?"
objectives:
- "Understand issues regarding string quoting, escaping and ..."
keypoints:
- "Add key points"
---

# Note - the Carpentries way seems to be that we intro quoting issues as they become relevant to the
# job in hand. So maybe do it that way? I think intro prameters here and quoting follows, then
# configuration from that.

Note - we probably want to introduce parameters first. The 'k' example in my assembly was perfect for this
but I can work out something else. Maybe a subsampling number? Or a quality cutoff?

Consider the following simple shell command:

  $ echo How "cool" are you?
  How cool are you?

The message is printed back, but not before the shell has interpreted the text as per the standard command-line
parsing rules

 * The quotes around "cool" have been removed
 * More subtly, the last work will be interpreted as a glob pattern...

  $ touch youX youY
  $ echo How "cool" are you?
  How cool are youX youY

Note: if you have certain shell settings you may see a warning about the unmatched glob pattern.

Anyone who has done any amount of programming or shell scripting will have come across quoting issues in code. Bugs related
to quoting can be particularly troublesome. In Snakemake these are particularly complex because the "shell" part of each rule
undergoes three rounds of interpretation before it is actually run:

 1. The string is parsed according to Python quoting rules
 1. Placeholders (eg. {input} {output}) are then replaced
 1. The resulting string goes to the shell and is subject to all shell parsing rules

In this episode we introduce the best practises for making your Snakefiles robust, and some simple rules to avoid most
mis-quoting complications.

Now have some exercises that illustrate the three points and introduce ways to deal with them:

r"""Strings like this""". Allows embedded newlines and literal \n \t or "quotes". Ie. Python does as little
as possible to the string. Good practise.

Double-escaping {{curlies}}. There is no good way to avoid this. Or is there? We can add them outside of the rule, maybe?
The AWK example is good here.

Shell parsing rules. Normally you want this because you want your pipes, redirects, shell variables, etc. to work.
If you choose filenames using only [0-9A-Z-a-z.\_] then you will always be "safe", but sometimes you have to deal with
awkwardly named files or parameters. If you want them unmolested. Use {foo:q}. Not '{foo}'!
Have an exercise that demonstrates why this is good.

How will I best demonstrate the power of {:q}?

letters = ['"', "'", "$ ", "\n\tx", "{:}"]
shell:
    "echo {letters:q}"

Hmmm. A little obscure. Well have a think.

{% include links.md %}

