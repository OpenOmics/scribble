<div align="center">
   
  <h1>scribble <sup>< / ></sup></h1>
  
  _**Scr**ipts/code-sn**i**ppets, somewhere in-**b**etween **b**ioinformatics and data science, and one**l**in**e**rs_
 
</div>

## Overview
Welcome to scribble! Scribble is a collection of bioinformatics and data science guides, scripts, code snippets, and oneliners. At the end of the day, we code to decipher the code of life. Did you write an awesome oneliner that became sentient? Scribble it down before it's lost forever!

### Hierarchy

Scribble is a collection of `guides`, `scripts`,  `snippets`, and `oneliners`. It is an amalgamation of different types and ideas. And as so, it is important to understand where you can find or save new information. Here is more information about each of the basic types that make up scribble. 

**Guides**

The guides directory of scribble contains stand-alone how-to guides. Guides can be as simple or elaborate as you want, but it is important to remember they are a set of human-readable instructions. The best guides leave nothing un-said, provide many examples, and assume the user knows nothing. 

**Scripts**

The scripts directory contains generalized scripts that accomplish a task, analysis, or objective. Did you write an awesome script that automates something cool? Did you create a script that produces a beautiful figure? Scribble it down here! If your script contains a few small input files (<200KB), please feel free to include them here as well. Scripts don't have to be written in python and R to be useful. Rmarkdown, quarto, jupyter notebooks, and even simple bash/sh scripts are all welcome! Please take some steps to ensure your scripts do not contain any hard-coded paths to files or resources. If your script has software requirements, please list them in the README accompanying your script.

**Snippets**

The snippets directory contains [code-snippets](https://en.wikipedia.org/wiki/Snippet_(programming)). These are usually small stand-alone pieces of reusable code. Snippets can act like templates for repeating certain acts/patterns/tasks. Did you figure out how to do something really cool in R or python that you want to share? Scribble it down here!

**One-liners**

The onliners directory contains a collection of oneliners. These are helpful shell commands-- that may or may not make blow your mind. Did you write a oneliner that automates your last job? Did your perl oneliner regex talk back to you? Scribble it down here!

## Contribute
To contribute to scribble, please create a new directory in one of the directories described above:
  
1. Add a README.md containing basic information about the script or code snippet
    - If your code relies on any third-party software or packages, please add a list of the required dependencies.  
2. Lint your scripts, code snippets, or config files and make changes as needed (optional, but recommended). Google maintains a [set of style guides](https://google.github.io/styleguide/) for different programming languages. These are the best practices developers/teams at Google follow. We don't enforce their recommendations; however, I would recommend reading through the style guides for python, R, and shell scripts at least once. There is a lot of great information in each of these guides!
    - `bash/sh`: [shellcheck](https://www.shellcheck.net/)
    - `python`: [flake8](https://github.com/pycqa/flake8), [pylint](https://pylint.pycqa.org/en/latest/), [black](https://black.vercel.app/)
    - `R`: [lintr](https://cran.r-project.org/web/packages/lintr/readme/README.html)
    - `JSON`: [jsonlint](https://jsonlint.com/)
    - `YAML`: [yamllint](http://www.yamllint.com/)
3. Upload any relevant files, scripts, or code snippets.
    - Please take some time to ensure your scripts or code snippets do ***not*** contain any hardcoded paths.
