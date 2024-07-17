---
title: Setup
---

## Software installation

These instructions set out how to obtain and install the course pre-requisites for Linux. It is
assumed that you have:

- access to the Bash shell on a fairly modern Linux system
- sufficient disk space (~1GB) to store the software and data

You **do not** need root/administrator level access.

:::::::::::::::::::::::::::::::::::::::::  callout

## Note

The recommended installation method here uses a *frozen* conda environment so that you will be
running the exact version of Snakemake, and other tools, that has been tested with this material.

There are other ways to install these tools, but Snakemake in particular is under active
development and so newer or older versions may not behave the same way that the material shows.

The frozen environment in `conda_env.yaml` only works on Linux just now, but we would welcome
[contributions](https://github.com/carpentries-incubator/snakemake-novice-bioinformatics/blob/gh-pages/CONTRIBUTING.md)
in this regard, such as setup instructions for Mac users. See the *Keeping the course updated*
section for more details.

::::::::::::::::::::::::::::::::::::::::::::::::::

We can install all the software packages, including Snakemake and the bioinformatics tools, via
Conda. For more info on Conda, see the first part of [episode 10
](episodes/10-conda_integration.md).

- If you don't already have Conda, start with [the Miniforge installer
  ](https://github.com/conda-forge/miniforge).

These are suggested commands to install and initialise Miniforge in a Linux Bash environment.

```bash
$ wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O installer.sh
$ bash ./installer.sh -b -u -p ~/miniforge3
$ rm installer.sh
$ ~/miniforge3/bin/conda init bash
$ exit
```

Then open a new shell. Your shell prompt should now start with `(base)`.

Now save the following into a file named `~/.condarc`:

```source
channels:
  - conda-forge
  - bioconda
solver: libmamba
channel_priority: strict
```

:::::::::::::::: callout

## Notes

- You may want to enable additional channels but all the packages needed in the course are in
  the two listed above.
- The use of *strict* channel priority when resolving dependencies is
  recommended by both Snakemake and [by Conda itself](
  https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html#strict).
- Setting `solver: libmamba` is nowadays preferred to explicitly running the `mamba` command.

:::::::::::::::::::::::::

Once this is set up, get [the *conda\_env.yaml* environment file](files/conda_env.yaml) and
create then activate the environment by running:

```bash
$ conda env update --file conda_env.yaml
$ conda activate snakemake_carpentry
```

You will also need a text editor such as [GEdit](https://help.gnome.org/users/gedit/stable/)
or [Nano](https://www.nano-editor.org/). These are standard on most Linux distros. Other text
editors will work fine but we only provide specific setup instructions for GEdit and Nano here.

## Obtaining the data

Download and unpack the sample dataset tarball from
[https://figshare.com/ndownloader/files/42467370](https://figshare.com/ndownloader/files/42467370)

You may do this in the shell with the command:

```bash
$ wget --content-disposition https://figshare.com/ndownloader/files/42467370
```

The [tar file](https://www.gnu.org/software/tar/manual/html_node/Tutorial.html)
needs to be unpacked to yield the directory of files used in the course. In the shell you may
do this with:

```bash
$ tar -xvaf data-for-snakemake-novice-bioinformatics.tar.xz
```

[See this link](https://figshare.com/articles/dataset/data-for-snakemake-novice-bioinformatics_tar_xz/19733338/1)
for details about this dataset and the redistribution license.

## Preparing your editor

There are some settings you should change in your editor to most effectively edit Snakemake
workflows. These are also good for editing most other types of script and code.

### GEdit

GEdit is the text editor that comes with the GNOME desktop and is a simple to use general-purpose
editor. Within the application menus it is normally just called "Text Editor" but you can also
start it from your shell by typing "`gedit`".

:::::::::::::::: callout

## Running gedit from the shell

If you type "`gedit`" in the shell, and GEdit is already running, a new tab will be created
in the existing editor and your shell prompt will return. If GEdit was not already running the
shell prompt in your terminal will not come back until GEdit is closed.

You can type "`gedit &`" to
get your shell prompt back immediately, but the program has a habit of printing sporadic warnings
into that terminal, so the cleanest option is to just start a new terminal window.

:::::::::::::::::::::::::

Snakemake uses Python file structure which is very fussy about the use of tab characters and line
breaks. Before writing any code in GEdit, you need to go into the preferences and select the
following settings:

- Insert spaces instead of tabs (under the "Editor" tab)
- Disable text wrapping (under the "View" tab)

The following settings are recommended but not required:

- Set the Syntax to "Python3", rather than "Plain Text"
- Set Tab width to 4
- Enable automatic indentation
- Highlight matching brackets
- Display line numbers

### Nano

The Nano editor works directly in the terminal and is found on virtually every Linux system.

The following command will start editing with the suggested settings, in particular regarding Tab
handling as mentioned above, and with Python syntax highlighting.

```bash
$ nano -wiSOE -T 4 -Y python Snakefile
```

To avoid typing all those options each time, you can add defaults to your `~/.nanorc` file:

```bash
$ nano ~/.nanorc
```

```source
set nowrap
set autoindent
set tabstospaces
set tabsize 4
set smooth
set morespace
```

Some of these are already defaults in later versions of Nano, but it doesn't hurt to have them
here anyway. There isn't a way to set the default syntax highlighting within this file.

