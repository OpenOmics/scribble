#!/usr/bin/env Rscript
# setup_renv.R
# Initialise a fresh renv in the current working directory, pin known-good
# package versions for the ShinyCell2 app, and install the remaining required
# packages at whatever version renv resolves.
#
# Run from the app's working directory:
#   Rscript setup_renv.R

# ---- 0. Bootstrap renv itself ------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  message(">> Installing renv from CRAN")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# ---- 1. Initialise a bare renv in the current directory ----------------------
# bare=TRUE means: create the structure but don't auto-discover/install anything.
# We want full control over what gets installed and at which version.
if (!file.exists("renv.lock")) {
  message(">> Initialising new renv in: ", normalizePath("."))
  renv::init(bare = TRUE, restart = FALSE)
} else {
  message(">> renv.lock already exists; using existing project")
  # Make sure renv is active in this session
  source("renv/activate.R")
}

# ---- 2. Pinned versions ------------------------------------------------------
# These are the must-pin and should-pin packages with specific versions known
# to be compatible with the ShinyCell2 codebase.
# ggplot2 + scales are the critical pair — ggplot2 4.x breaks the +.gg chains
# used throughout shinyFunc.R via the S7 migration.
pinned <- c(
  "ggplot2@3.5.2",      # MUST be 3.x — 4.0 breaks +.gg dispatch
  "scales@1.3.0",       # matched to ggplot2 3.5.x
  "S7@0.2.0",           # belt-and-suspenders against ggplot2 4 sneaking in
  "ggrepel@0.9.5",      # geom_text_repel(bg.color=, bg.r=) args
  "ggpubr@0.6.0",       # stat_compare_means signature
  "ggdendro@0.2.0",     # dendro_data() return shape used in scBubbHeat
  "hdf5r@1.3.11",       # H5File$new(..., mode="r") + $read(args=list(...))
  "data.table@1.16.4",  # heavy [.data.table use, :=, dcast.data.table, CJ
  "Matrix@1.6-5",       # sparse matrix infrastructure; pin major version
  "DT@0.33",            # datatable() extensions="Buttons" + formatRound
  "shiny@1.8.1.1",      # shinyServer/shinyUI, brushOpts, downloadHandler
  "shinyhelper@0.3.2",  # observe_helpers() + helper() pipe-chained
  "magrittr@2.0.3"      # %>%
)

# ---- 3. Unpinned but required ------------------------------------------------
# Resolved to whatever renv pulls in (latest compatible version).
# gridExtra: grid.arrange/arrangeGrob in scBubbHeat
unpinned <- c(
  "gridExtra"
)

# ---- 4. Install the pinned versions first ------------------------------------
# Doing pinned first means later unpinned installs can't pull in a newer
# conflicting version of something we care about — renv resolves against
# what's already installed.
message("\n>> Installing pinned packages")
for (pkg in pinned) {
  message(sprintf("   - %s", pkg))
  tryCatch(
    renv::install(pkg, prompt = FALSE),
    error = function(e) {
      stop(sprintf("Failed to install '%s': %s", pkg, conditionMessage(e)))
    }
  )
}

# ---- 5. Install the remaining required packages ------------------------------
message("\n>> Installing unpinned required packages")
for (pkg in unpinned) {
  pkg_name <- sub("@.*$", "", pkg)
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    message(sprintf("   - %s (already installed)", pkg_name))
  } else {
    message(sprintf("   - %s", pkg))
    tryCatch(
      renv::install(pkg, prompt = FALSE),
      error = function(e) {
        stop(sprintf("Failed to install '%s': %s", pkg, conditionMessage(e)))
      }
    )
  }
}

# ---- 6. Snapshot to the lockfile ---------------------------------------------
# type="all" captures every installed package, even those not directly
# referenced from source files. Use type="implicit" instead if you want renv
# to scan code and only lock what it can find references to.
message("\n>> Writing renv.lock")
renv::snapshot(prompt = FALSE, type = "all")

# ---- 7. Verify -------------------------------------------------------------
message("\n>> Verification")
required <- c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr",
              "ggplot2", "ggrepel", "hdf5r", "ggdendro", "gridExtra", "ggpubr",
              "scales", "S7")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("The following packages failed to install: ",
       paste(missing, collapse = ", "))
}

# Report installed versions for the record
cat("\nInstalled versions:\n")
for (pkg in required) {
  cat(sprintf("  %-15s %s\n", pkg, as.character(packageVersion(pkg))))
}

# ---- 8. Smoke test for the ggplot2 issue -------------------------------------
# Confirm the specific +.gg failure mode is resolved before declaring success.
message("\n>> ggplot2 smoke test")
suppressPackageStartupMessages(library(ggplot2))
p <- tryCatch(
  ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_bw(),
  error = function(e) {
    stop("ggplot2 smoke test FAILED: ", conditionMessage(e),
         "\n   The renv has the same problem as before. Check `packageVersion('ggplot2')`.")
  }
)
message("   ggplot2 smoke test passed")

message("\nDone. The app should now run cleanly with:  shiny::runApp()")