---
title: Setup
---

These instructions set out how to obtain and install the course pre-requisites for Linux. It is assumed
that you have:
* access to the Bash shell on a fairly modern Linux system
* sufficient disk space (~1GB) to store the software and data

You **do not** need root/administrator level access.

> ## Note
> There are other ways to install these tools, and it should also be possible to run the exercises on
> non-Linux systems. However, this is not covered here.
{: .callout}

## Software and data installation

To prepare for the course:

* Download and unpack [the sample dataset tarball](https://ndownloader.figshare.com/files/28531743)
* Install software packages, including Snakemake and the bioinformatics tools, via Conda:
  * If you don't already have Conda, get [the Miniconda installer](https://docs.conda.io/en/latest/miniconda.html) and follow instructions
  * Get [the environment file](files/conda_env_all.yaml) and run `conda env update --file conda_env_all.yaml`
  * Ensure the right environment is active - `conda activate snakemake_dash`
* Set up the GEdit or Nano text editor (other text editors will work fine but we only provide specific setup instructions
  for GEdit and Nano).

## Preparing your editor

### GEdit

GEdit is the text editor that comes with the GNOME desktop and is a great general-purpose editor. Within the application
menus it is normally just called "Text Editor" but you can also start it from the shell by typing "gedit &".

Snakemake uses Python file structure which is very fussy about editor settings, particularly the use of Tab characters.
Before using GEdit, you need to go into the preferences and select the following settings:

* Insert spaces instead of tabs (under the "Editor" tab)
* Disable text wrapping (under the "View" tab)

The following settings are also recommended:

* Set the Syntax to "Python3", rather than "Plain Text"
* Set Tab width to 4
* Enable automatic indentation
* Hilight matching brackets
* Display line numbers

### Nano

The Nano editor works directly in the terminal and is found on virtually every Linux system.

The following command will start editing with the suggested settings, in particular regarding Tab handling as mentioned
above.

```
$ nano -wiSOE -T 4 -Y python Snakefile
```

You probably want to make an alias for this:

```
$ alias sfedit="nano -wiSOE -T 4 -Y python"
```

{% include links.md %}
