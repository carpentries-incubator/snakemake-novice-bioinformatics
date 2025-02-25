---
title: Processing lists of inputs
teaching: 50
exercises: 30
---

::::::::::::::::::::::::::::::::::::::: objectives

- Use Snakemake to process all our samples at once
- Make a summary of all read counts

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- How do I define a default set of outputs for my Snakefile?
- How do I make rules which combine whole lists of files?
- How do I process all available input files at once?

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is the Snakefile](files/ep04.Snakefile) you should have to start the
episode. We didn't modify it during the last episode.*

## Defining a list of samples to process

So far, we've told Snakemake what output files to generate by giving the names of the desired files
on the command line. Often you want Snakemake to process all the available samples. How can we do
this?

The `yeast/reads` directory contains results from three conditions: `ref`, `etoh60` and `temp33`.
There are three replicates for each condition. There is minor inconsistency in the naming
convention, as the `ref` files do not have an underscore before the replicate number. Consistent
naming is important for Snakemake, so let's fix up the names before we go any further.

A good way to do this is by making symlinks, because that way you don't lose sight of the original
file names.

```bash
$ mv reads original_reads
$ mkdir reads
$ cd reads
$ ln -s ../original_reads/* .
$ rename -v -s ref ref_ *
$ cd ..
```

:::::::::::::::::::::::::::::::::::::::::  callout

## File renaming

The **rename** command here is [the one provided by Bioconda
](https://plasmasturm.org/code/rename/). Other Linux systems may have a different rename
command installed by default.

::::::::::::::::::::::::::::::::::::::::::::::::::

Having harmonized the file names we'll tell Snakemake about our conditions
and replicates. To do this, we can define some lists as Snakemake **global variables**.

Global variables should be added before the rules in the Snakefile.

```source
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]
```

- Unlike with variables in shell scripts, we can put spaces around the `=` sign, but they are
  not mandatory.
- The lists of quoted strings are enclosed in square brackets and comma-separated. If you know any
  Python you'll recognise this as Python list syntax.
- A good convention is to use capitalized names for these variables, but this is not mandatory.
- Although these are referred to as variables, you can't actually change the values once the
  workflow is running, so lists defined this way are more like constants.

## Using a Snakemake rule to define a batch of outputs

We'll add another rule to our Snakefile. This special target rule will have an `input` section but
no `output` or `shell` sections (yet).

```source
rule all_differences:
    input: expand("{cond}_{rep}_1.reads_removed.txt", cond=CONDITIONS, rep=REPLICATES)
```

The `expand(...)` function in this rule generates a list of filenames, by taking the first thing in
the parentheses as a template and replacing `{cond}` with all the `CONDITIONS` and `{rep}` with all
the `REPLICATES`. Since there are 3 of each, this will yield 9 combinations - ie. 9 files we want
to make.

This list goes into the *input* section of the rule. You might think that since these filenames are
outputs that we want from our workflow they should go into the *output* section. However, remember
that *output*s of a rule are things the rule can make itself, and this rule doesn't actually make
anything. It's just a placeholder for a bunch of filenames.

We now tell Snakemake to make all these files by using the *target rule name* on the command line:

```bash
$ snakemake -j1 -p all_differences
```

Here, Snakemake sees that *all\_differences* is the name of a rule in the Snakefile, so rather than
trying to make a file literally named `all_differences` it looks at all the input files for the
rule and tries to make them. In this case, all of the inputs to *all\_differences* can be made by
the *calculate\_difference* rule, and all of the inputs for those jobs are made by *trimreads* and
*countreads*. The resulting workflow is the same as if we had typed out all 9 of the filenames on
the command line.

If you don't specify a target rule name or any file names on the command line when running
Snakemake, the default is to use **the first rule** in the Snakefile as the target. So if
*all\_differences* is defined at the top, before the other rules, you can simply say:

```bash
$ snakemake -j1 -p
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Rules as targets

Giving the name of a rule to Snakemake on the command line only works when that rule has
*no wildcards* in the outputs, because Snakemake has no way to know what the desired wildcards
might be. You will see the error "Target rules may not contain wildcards." This can also happen
when you don't supply any explicit targets on the command line at all, and Snakemake tries to run
the first rule defined in the Snakefile.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## Counting all the reads

Check that the *all\_differences* rule is working. Now adapt the Snakefile so that it makes all the
*reads_removed.txt* for both of the pairs of reads (`_1.fq` and `_2.fq`). So you should end up
with 18 output files in total.

:::::::::::::::  solution

## Solution

This will work.

```source
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]
READ_ENDS  = ["1", "2"]

# Rule to make all reads_removed.txt files at once
rule all_differences:
    input: expand("{cond}_{rep}_{end}.reads_removed.txt", cond=CONDITIONS, rep=REPLICATES, end=READ_ENDS)
```

Alternatively you can put the list directly into the `expand()` function rather than declaring
more variables. To aid readability of the code it's also possible to split the function over
more than one line, but note that this only works if you put a newline after the `input:`
line too.

```source
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

# Rule to make all reads_removed.txt files at once
rule all_differences:
    input:
        expand( "{cond}_{rep}_{end}.reads_removed.txt", cond = CONDITIONS,
                                                        rep  = REPLICATES,
                                                        end  = READ_ENDS )
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Rules that combine multiple inputs

Our *all\_differences* rule is a rule which takes a list of input files. The length of that list is
not fixed by the rule, but can change based on `CONDITIONS` and `REPLICATES`. If we want to perform
some combining operation on the whole list of files, we can add `output` and `shell` sections to
this rule.

In typical bioinformatics workflows, the final steps will combine all the results together into
some big report. For our final workflow we'll be doing this with *MultiQC*, but as a simple first
example, let's just concatenate all the text files. In the shell this would be:

```bash
$ cat *.reads_removed.txt > all_reads_removed.txt
```

In the Snakemake rule we just say:

```bash
shell:
  "cat {input} > {output}"
```

Within a rule definition you can combine named inputs and list inputs - any named input can be list
of files rather than just a single file. When you use the `{input.name}` placeholder in the shell
command it will expand to the full list.

:::::::::::::::::::::::::::::::::::::::  challenge

## Combining the inputs of the *all\_differences* rule

1. Make it so that the *all\_differences* rule concatenates all the `.reads_removed.txt` files into
   a single output file named `all_reads_removed.txt`.
2. Adapt the rule further so that there are two outputs named `all_read1_removed.txt`
   and `all_read2_removed.txt`, and the respective counts go into each.

*Hint: you can put two commands into the `shell` part, separated by a semicolon `;`.*

:::::::::::::::: solution

## Solution

```source
rule all_differences:
    input:
        expand( "{cond}_{rep}_{end}.reads_removed.txt", cond = CONDITIONS,
                                                        rep  = REPLICATES,
                                                        end  = READ_ENDS )
    output:
        "all_reads_removed.txt"
    shell:
        "cat {input} > {output}"
```

And for the version with separate output files.

```source
rule all_differences:
    input:
        read1 = expand( "{cond}_{rep}_1.reads_removed.txt", cond  = CONDITIONS,
                                                            rep   = REPLICATES ),
        read2 = expand( "{cond}_{rep}_2.reads_removed.txt", cond  = CONDITIONS,
                                                            rep   = REPLICATES ),
    output:
        read1 = "all_read1_removed.txt",
        read2 = "all_read2_removed.txt",
    shell:
        "cat {input.read1} > {output.read1} ; cat {input.read2} > {output.read2}"
```

This works, but may more elegantly be broken down into two rules, so that once again the target
rule has only *input* and no *outputs*. There is a quirk where we
need to put extra braces around `{{end}}` in the place where this is a wildcard, and not to be
expanded out.

```source
rule all_differences:
    input: expand("all_read{end}_removed.txt", end=["1","2"])

rule all_differences_per_end:
    input:
        expand( "{cond}_{rep}_{{end}}.reads_removed.txt", cond  = CONDITIONS,
                                                          rep   = REPLICATES )
    output:
        "all_read{end}_removed.txt"
    shell:
        "cat {input} > {output}"
```

To run any version of the rule:

```bash
$ snakemake -j1 -p all_differences
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

## Dynamically determining the inputs

In the shell we can match multiple filenames with **glob patterns** like `original_reads/*` or
`reads/ref?_?.fq`. Snakemake allows you to do something similar with the `glob_wildcards()`
function, so we'll use this in our Snakefile.

```source
CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition
print("Conditions are: ", CONDITIONS)
```

Here, the list of conditions is captured from the files seen in the reads directory. The pattern
used in `glob_wildcards()` looks much like the input and output parts of rules, with a wildcard in
{curly brackets}, but here the pattern is being used to search for matching files. We're only
looking for read 1 of replicate 1 so this will return just three matches

Rather than getting the full file names, the function yields the values of the wildcard, which we
can assign directly to a list. The `print()` statement will print out the value of `CONDITIONS`
when the Snakefile is run (including dry-run mode, activated with the `-n` option as mentioned in
[Episode 2](02-placeholders.html)), and reassures us that the list really is the same as before.

```bash
$ snakemake -j1 -F -n -p all_differences
Conditions are:  ['etoh60', 'temp33', 'ref']
Building DAG of jobs...
Job stats:
job                        count
-----------------------  -------
all_differences                1
all_differences_per_end        2
calculate_difference          18
countreads                    36
trimreads                     18
total                         75
...
```

Using `glob_wildcards()` gets a little more tricky when you need a more complex match. To refine
the match we can quickly test out results by activating the Python interpreter. This saves editing
the Snakefile and running Snakemake just to see what `glob_wildcards()` will match.

```python
$ python3
>>> from snakemake.io import glob_wildcards, expand
>>> glob_wildcards("reads/{condition}_1_1.fq")
Wildcards(condition=['etoh60', 'temp33', 'ref'])
```

This is the result we got before. So far, so good.

:::::::::::::::::::::::::::::::::::::::::  callout

## The Python interpreter

The Python interpreter is like a special shell for Python commands, and is a familiar to anyone
who has learned Python. Because Snakemake functions are actually Python functions they can be run
from the interpreter after being imported.

Note that `>>>` is the Python interpreter prompt, and not part of the command to be typed. You
can exit back to the regular shell prompt by typing `exit()`. If you do exit and then restart the
interpreter, you will need to repeat the `import` line.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

## 'Globbing' the list of samples

Staying in the Python interpreter, use `glob_wildcards()` to list the names of all nine samples,
that is the three replicates of each condition.

:::::::::::::::  solution

## Solution

```python
>>> glob_wildcards("reads/{sample}_1.fq")
Wildcards(sample=['temp33_3', 'temp33_2', 'etoh60_1', 'etoh60_3', 'ref_2', 'temp33_1', 'etoh60_2', 'ref_1', 'ref_3'])
```

:::::::::::::::::::::::::



Still in the Python interpreter, use the `expand()` function in combination with the
`glob_wildcards()` to make a list of all the count files, and then all the *kallisto\_quant*
output files that would be generated for all the nine samples.

Remember you can save the result of `glob_wildcards()` to a variable - this works the same way in
the Python interpreter as it does in the Snakefile.

:::::::::::::::  solution

## Solution

```python
>>> SAMPLES = glob_wildcards("reads/{sample}_1.fq").sample
>>> expand("{indir}.{sample}_{end}.fq.count", indir=['reads', 'trimmed'], sample=SAMPLES, end=["1","2"])
...
>>> expand("kallisto.{sample}/{outfile}", sample=SAMPLES, outfile=['abundance.h5', 'abundance.tsv', run_info.json'])
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Glob with multiple wildcards

If there are two or more wildcards in the glob pattern, dealing with the result becomes a little
more tricky, but it does work.

Here is one way to re-combine two wildcards using `expand()` and `zip` (demonstrated in the
Python interpreter as above).

```python
$ python3
>>> from snakemake.io import glob_wildcards, expand
>>> DOUBLE_MATCH = glob_wildcards( "reads/{condition}_{samplenum}_1.fq" )
>>> DOUBLE_MATCH
Wildcards(condition=['temp33', 'temp33', 'etoh60', 'etoh60', 'ref', 'temp33', 'etoh60', 'ref', 'ref'],
          samplenum=['3', '2', '1', '3', '2', '1', '2', '1', '3'])
>>> SAMPLES = expand(
...     "{condition}_{samplenum}",
...     zip,
...     condition = DOUBLE_MATCH.condition,
...     samplenum = DOUBLE_MATCH.samplenum )
>>> SAMPLES
['temp33_3', 'temp33_2', 'etoh60_1', 'etoh60_3', 'ref_2', 'temp33_1', 'etoh60_2', 'ref_1', 'ref_3']
>>> sorted(SAMPLES)
['etoh60_1', 'etoh60_2', 'etoh60_3', 'ref_1', 'ref_2', 'ref_3', 'temp33_1', 'temp33_2', 'temp33_3']
```

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

## Rules that make multiple outputs

If we can have rules that combine lists of files, can we do the opposite and have a rule that
produces a list of outputs?

The answer is yes, but the situation is not completely symmetrical. Remember that Snakemake works
out the full list of input and output files to every job *before* it runs a single job in the
workflow. For a combining step, Snakemake will know how many samples/replicates are being
combined from the start. For a splitting step, it may or may not be possible to predict the
number of output files in advance. For cases where you really do need to handle a dynamic list of
outputs, Snakemake has things called
[checkpoint rules](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution).

In practise these are very rarely needed, so we'll not be covering them here in the course.

::::::::::::::::::::::::::::::::::::::::::::::::::

*For reference, [this is a Snakefile](files/ep06.Snakefile) incorporating the changes made in
this episode.*



:::::::::::::::::::::::::::::::::::::::: keypoints

- Rename your input files if necessary to maintain consistent naming
- List the things you want to proces as global variables, or discover input files with
  `glob_wildcards()`
- Use the `expand()` function to generate lists of filenames you want to combine
- These functions can be tested in the Python interpreter
- Any `{input}` to a rule can be a variable-length list
- (But variable lists of outputs are trickier and rarely needed)

::::::::::::::::::::::::::::::::::::::::::::::::::

