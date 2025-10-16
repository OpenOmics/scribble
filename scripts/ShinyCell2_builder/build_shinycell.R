#!/usr/bin/env Rscript

# Misc helper functions
options(error = function() traceback(2))
err <- function(...) {
  cat(sprintf(...), sep = "\n", file = stderr())
}
fatal <- function(...) {
  err(...)
  quit(status = 1)
}

# get and check libraries
required_packages <- c("Seurat", "optparse", "remotes")
required_gh_packages <- c(ShinyCell2 = "OpenOmics/ShinyCell2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in names(required_gh_packages)) {
  if (!require(pkg, character.only = TRUE)) {
    remotes::install_github(required_gh_packages[[pkg]])
  }
}

# Load libraries
suppressPackageStartupMessages(library(ShinyCell2, quietly = TRUE))
suppressPackageStartupMessages(library(Seurat, quietly = TRUE))
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

# Command line arguments
option_list <- list(
  make_option(
    c("-j", "--obj"),
    dest = "object",
    metavar = "Seurat Object file path",
    type = "character",
    default = NULL,
    help = "RDS file created with saveRDS containing a seurat object"
  ),
  make_option(
    "--proj",
    dest = "project",
    type = "character",
    metavar = "Project Name",
    default = NULL,
    help = "Project name, becomes title of the app, example: NCBR-34"
  ),
  make_option(
    "--rmmeta",
    dest = "meta.to.rm",
    metavar = "comma,delimited,metadata,column,names",
    type = "character",
    default = NULL,
    help = paste0(
      "Comma delimited list of names in seurat_object@meta.data ",
      "that will be dropped prior to deploying to ShinyCell2"
    )
  ),
  make_option(
    "--defred",
    type = "character",
    metavar = "seurat_reduction_name",
    dest = "defaultreduction",
    default = NULL,
    help = paste0(
      "The default reduction to use with ShinyCell2, ",
      "this value must exist in Seurat::DefaultDimReduc(obj).\n The ",
      "two major prinipal components for this reductions must ",
      "be labeled the same label with a 1 and a 2 trailing it.\n ",
      "i.e. defaultreduction = UMAP, UMAP1 and UMAP2 are the two components ",
      "that must exist"
    )
  ),
  make_option(
    c("-l", "--maxlevels"),
    type = "integer",
    dest = "max.levels",
    metavar = "MAX LEVELS [int]",
    default = NULL,
    help = paste(
      "The maximum allowable amount levels/factors of categorical values,",
      "any seurat meta.data with factors/levels greater than",
      "this number will be dropped. Default is 50. ",
      "i.e. if seurat_object@meta.data$seurat_clusters has clusters",
      "1->120 (120 factors/levels), and max levels is 50 this metadata",
      "will be discarded.",
      sep = " "
    )
  ),
  make_option(
    c("-a", "--assay"),
    type = "character",
    dest = "assaytouse",
    metavar = "ASSAY_NAME [str]",
    default = NULL,
    help = paste(
      "The assay to utilize for ShinyCell2 web application.",
      "Comma delimit multiple assays in a single string e.g.: RNA,spatial,ATAC,etc.",
      "This will default to the first assay in the Seurat object (object@assays)",
      sep = " "
    )
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

# setup opt parse variables for downstream usage into shinycell2
rds_file <- opt$object
project_name <- opt$project
required_args <- c("object", "project")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]
shiny_app_dir <- file.path(paste0(getwd(), '/shinyApp'))

if (length(missing_args) > 0) {
  cat("Error: Missing required arguments:", paste(missing_args, collapse = ", "), "\n\n")
  print_help(opt)
  quit(status = 1)
}

seurat_obj <- readRDS(rds_file)
if (is.null(opt$meta.to.rm) || is.na(opt$max.levels) || opt$meta.to.rm == "" || opt$meta.to.rm == "NA") {
  rmmeta <- NULL
} else {
  if ("," %in% opt$meta.to.rm) {
    rmmeta <- unlist(strsplit(opt$meta.to.rm, ",", fixed = TRUE))
  } else {
    rmmeta <- c(trimws(gsub("[\r\n]", "", opt$meta.to.rm)))
  }
}
if (is.null(opt$assaytouse) || is.na(opt$assaytouse) || opt$assaytouse == "" || opt$assaytouse == "NA") {
  assaytouse <- NULL
} else {
  if ("," %in% opt$assaytouse) {
    assaytouse <- unlist(strsplit(opt$assaytouse, ",", fixed = TRUE))
  } else {
    assaytouse <- c(trimws(gsub("[\r\n]", "", opt$assaytouse)))
  }
}
if (is.null(opt$defaultreduction) || is.na(opt$defaultreduction) || opt$defaultreduction == "" || opt$defaultreduction == "NA") {
  defaultreduction <- NULL
} else {
  if (opt$defaultreduction %in% names(seurat_obj@reductions)) {
    this_key <- seurat_obj@reductions[[opt$defaultreduction]]@key
    defaultreduction <- c(paste0(this_key, "1"), paste0(this_key, "2"))
  } else {
    fatal(paste0("`", opt$defaultreduction, "` reduction not found in seurat object!"))
  }
}
if (!is.null(opt$max.levels) || !is.na(opt$max.levels) || opt$max.level == "" || opt$max.level == "NA") {
  max.levels <- opt$max.levels
} else {
  max.levels <- NULL
}
dir.create(shiny_app_dir, showWarnings = FALSE)

# Sanity check: Does the RDS file
# actually contain a seurat object?
if (class(seurat_obj) == "SeuratObject") {
  # It doesn't look like it...
  err("Fatal Error: Failed to provide an RDS file with a Seuart Object.")
  fatal(" └── Please create a new RDS file with a seurat object!")
}

# Remove unsupported assay
# ShinyCell2 supports:
#   - CITEseq
#   - spatial
#   - scATAC (Signac)
# ShinyCell2 does NOT support:
#   - celltype "assays" added by azimuth
#   - HTO
for (assay in Assays(seurat_obj)) {
  if (startsWith(assay, "prediction.score.celltype.l")) {
    seurat_obj[[assay]] <- NULL
  }
}

unsupported_assays <- c("HTO")
for (assay in unsupported_assays) {
  if (assay %in% names(seurat_obj@assays)) {
    cat(paste0(assay, "unsupported assay removed!", sep = " "))
    seurat_obj[[assay]] <- NULL
  }
}

# Create ShinyCell config file
# to make the application
config_params <- list()
if (!is.null(max.levels)) {
  config_params$maxLevels <- max.levels
}

shinycell_config <- do.call(
  createConfig,
  c(seurat_obj, config_params)
)

remove_metas <- c()

if (!is.null(rmmeta)) {
  remove_metas <- c(remove_metas, rmmeta)
}

for (config_label in shinycell_config$ID) {
  for (assay in unsupported_assays) {
    if (grepl(assay, config_label, fixed = TRUE)) {
      remove_metas <- c(remove_metas, config_label)
    }
  }
}

shinycell_config <- delMeta(shinycell_config, remove_metas)

# Build the Shiny Application,
# in the default location for
# Shiny/Posit server: i.e.
# /srv/shiny-server/${app_name}
files_params <- list(
  seurat_obj,
  shinycell_config,
  shiny.dir = shiny_app_dir,
  shiny.prefix = "sc1"
)

if (!is.null(defaultreduction)) {
  files_params$dimred.to.use <- opt$defaultreduction
  files_params$default.dimred <- defaultreduction
}

do.call(
  makeShinyFiles,
  files_params
)
makeShinyCodes(
  shiny.title = project_name,
  shiny.dir = shiny_app_dir,
  shiny.prefix = "sc1"
)

