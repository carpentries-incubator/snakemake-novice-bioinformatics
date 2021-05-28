# Snakemake for Bioinformatics

TODO - get this listed on the  [Community Developed Lessons page][community-lessons] once ready.

## Snakemake for Bioinformatics

This lesson introduces the Snakemake workflow system in the context of a bioinformatics data
analysis task.

To quote from the [official Snakemake documentation](https://snakemake.readthedocs.io/):

> The Snakemake workflow management system is a tool to create reproducible and scalable data analyses.
> Workflows are described via a human readable, Python based language. They can be seamlessly scaled to
> server, cluster, grid and cloud environments, without the need to modify the workflow definition.
> Finally, Snakemake workflows can entail a description of required software, which will be automatically
> deployed to any execution environment.

Snakemake originated as, and remains most popular as, a tool for bioinformatics. This is how we present
it here. However, Snakemake is a general-purpose system and may be used for all manner of data processing
tasks.

Snakemake is a superset of the Python language and as such can draw on the full power of Python, but you
do not need to be a Python programmer to use it. This leson **assumes no prior knowledge of Python** and
intruduces just a few concepts as needed to construct useful workflows.

## Lesson development checklist

Before you begin developing your new lesson,
here are a few things we recommend you do:

* [✓] Decide on a title for your new lesson!
  Once you've chosen a new title, you can set the value for `lesson_title`
  in [`_config.yml`](_config.yml)
* [ ] Add the URL to your built lesson pages to the repository description\*
* [ ] [Add relevant topic tags to your lesson repository][cdh-topic-tags].
* [ ] Fill in the fields marked `FIXME` in:
  * this README
  * [`_config.yml`](_config.yml)
* [✓] If you're going to be developing lesson material for the first time
  according to our design principles,
  consider reading the [Carpentries Curriculum Development Handbook][cdh]
* [✓] Consult the [Lesson Example][lesson-example] website to find out more about
  working with the lesson template
* [✓] If you are planning to write your lesson in RMarkdown,
  [create a `gh-pages` branch and set this as the default branch in your repository settings][change-default-branch]
* [ ] Update this README with relevant information about your lesson
  and delete this section


\* To set the URL on GitHub, click the gear wheel button next to **About**
on the right of the repository landing page.
The lesson URL structure is **https://carpentries-incubator.github.io/<repository-slug\>**:
a repository at https://github.com/carpentries-incubator/new-lesson/ will have pages at
the lesson URL https://carpentries-incubator.github.io/new-lesson/.


## Contributing

We welcome all contributions to improve the lesson! Maintainers will do their best to help you if you have any
questions, concerns, or experience any difficulties along the way.

We'd like to ask you to familiarize yourself with our [Contribution Guide](CONTRIBUTING.md) and have a look at
the [more detailed guidelines][lesson-example] on proper formatting, ways to render the lesson locally, and even
how to write new episodes.

Please see the current list of [issues][FIXME] for ideas for contributing to this
repository. For making your contribution, we use the GitHub flow, which is
nicely explained in the chapter [Contributing to a Project](http://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) in Pro Git
by Scott Chacon.
Look for the tag ![good_first_issue](https://img.shields.io/badge/-good%20first%20issue-gold.svg). This indicates that the maintainers
will welcome a pull request fixing this issue.


## Maintainer(s)

Current maintainers of this lesson are

* [Tim Booth](http://github.com/tbooth)
* [Graeme Grimes](http://github.com/ggrimes)

## Authors

A list of contributors to the lesson can be found in [AUTHORS](AUTHORS)

## Citation

To cite this lesson, please consult with [CITATION](CITATION)

[cdh]: https://cdh.carpentries.org
[cdh-topic-tags]: https://cdh.carpentries.org/the-carpentries-incubator.html#topic-tags
[change-default-branch]: https://docs.github.com/en/github/administering-a-repository/changing-the-default-branch
[community-lessons]: https://carpentries.org/community-lessons
[lesson-example]: https://carpentries.github.io/lesson-example
