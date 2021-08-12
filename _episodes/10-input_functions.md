I think we get to this in regard to Salmon and Kallisto. Let me see...

OK, so if ref_1, ref_2 and ref_3 are really three replicates then we should be running them one at a time.
But if they are all reads from the same sample, we want to add them all at once.

So for this we'd run kallisto three times, and each time would have 9 inputs.

Like:

rule kallisto_quant:
    output:
        outdir = directory("kallisto.{sample}"),
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = [ "trimmed/{sample}_1_1.fq", "trimmed/{sample}_2_1.fq", "trimmed/{sample}_3_1.fq" ],
        fq2   = [ "trimmed/{sample}_1_2.fq", "trimmed/{sample}_2_2.fq", "trimmed/{sample}_3_2.fq" ],
    shell:
        "kallisto quant -i {input.index} -o {output.outdir} {input.fq1} {input.fq2}"

Hmmm. So because I always have 3 repeats I can get away without an input funtion.
I could use glob_wildcards() though and it would be legit.

Then how do I get the fq1 and fq2 into the shell command? It's tricky - I actually need to turn the shell
part into a run part:

run:
    foo = [ i for p in zip(input.list1, input.list2) for i in p ]
    shell("echo {foo}")

Oh gods we've hit advanced Python syntax. Is there a small-brain way to do this??? I think we just have to make
the fqs a single list in the first place, then use expand() in the input function. Yep. Do this. It computes.
