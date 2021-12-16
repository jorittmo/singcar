# Contributing to `singcar`

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

First of all, thanks for considering contributing to `singcar`!

`singcar` is an open source project, maintained by people who care. We are not directly funded to do so.

[repo]: https://github.com/jorittmo/singcar
[docs]: https://cran.r-project.org/web/packages/singcar/singcar.pdf
[issues]: https://github.com/jorittmo/singcar/issues
[new_issue]: https://github.com/jorittmo/singcar/issues/new
[email]: mailto:j.rittmo@gmail.com

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want to know more about why and how to contribute to open source projects like this one, see this [Open Source Guide](https://opensource.guide/how-to-contribute/).

### Share the love ❤️

Think `singcar` is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using `singcar` for a paper you are writing? Consider citing it, e.g. by using `utils::citation("singcar")`.

### Ask a question️

Using `singcar` and got stuck? Browse the [documentation][docs] or [vignette](https://cran.r-project.org/web/packages/singcar/vignettes/singcar_vignette.html) to see if you can find a solution. Still stuck? Post your question as an [issue on GitHub][new_issue]. I will try my best to address it. 

Want to ask a question in private? Contact the package maintainer by [email][email].

### Feature suggestions

Have an idea for a new `singcar` feature? Take a look at the [documentation][docs] and [issue list][issues] to see if it isn't included or suggested yet. If not, suggest your idea as an [issue on GitHub][new_issue].

### Report a bug or issue

Report bugs or issues as a new [new issue][new_issue] in the issue tracker, so we can fix it. A good bug report makes it easier for us to do so, so please include:

* Your session info (e.g. using `utils::sessionInfo()`).
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed (but minimal) steps to reproduce the bug, preferably with a reproducible example using [reprex](https://reprex.tidyverse.org/). 

### Improve the documentation

Noticed a typo in the vignette? Think a function could use a better example? Good documentation makes all the difference, so your help to improve it is very welcome!

#### Function documentation

Functions are described as comments near their code and translated to documentation using [`roxygen2`](https://klutometis.github.io/roxygen/). If you want to improve a function description:

1. Go to `R/` directory in the [code repository][repo].
2. Look for the file with the name of the function.
3. [Propose a file change](https://help.github.com/articles/editing-files-in-another-user-s-repository/) to update the function documentation in the roxygen comments (starting with `#'`).

### Contribute code

Care to fix bugs or implement new functionality for `singcar`? Have a look at the [issue tracker][issues], leave a comment on things you are
working on and see the pull request guidelines below.

# Pull request guidelines

We try to follow the [GitHub flow](https://guides.github.com/introduction/flow/) for development.

1. Fork [this repo][repo] and clone it to your computer. To learn more about this process, see [this guide](https://guides.github.com/activities/forking/).
2. If you have forked and cloned the project before and it has been a while since you worked on it, [pull changes from the original repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/) to your clone by using `git pull upstream master`.
3. Open the RStudio project file (`.Rproj`).
4. Make your changes:
    * Write your code.
    * Test your code (bonus points for adding unit tests).
    * Document your code (see function documentation above).
    * Check your code with `devtools::check()` and aim for 0 errors and warnings.
5. Commit and push your changes.
6. Submit a [pull request](https://guides.github.com/activities/forking/#making-a-pull-request).

### Reference

These contribution guidelines have been adopted and adapted from [here](https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c).
