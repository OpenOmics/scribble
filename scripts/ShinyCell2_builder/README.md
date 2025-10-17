# Introduction

This script is used to create the ShinyCell2 application for taring and usage on a local system

# Setup

To setup an environment to run this with I recommend using renv like this:

```R
install.packages('renv')
setwd('/path/to/where/seurat/object/is')
renv::init()
renv::load()
renv::install(c("Seurat", "optparse", "remotes"))
remotes::install_github("OpenOmics/ShinyCell2")
```

> [!IMPORTANT]
> This script is intended to be used with OpenOmics/ShinyCell2 fork, athough utility with the original shinycell2 is possible
