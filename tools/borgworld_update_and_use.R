# =========================================================
# borgworld: update package and use from client projects
# =========================================================

# ---- A0. Paths (edit these) ----
# Path to your borgworld package source
borgworld_path <- "C:/Users/sborg2/R/borgworld"

# ---- A1. Add or update functions (in borgworld project) ----
# Put these in borgworld/R/bsum.R (or similar), then save the file.
#' @export
bsum <- function(x, ...) sum(x, na.rm = TRUE, ...)
#' @export
bmean <- function(x, ...) mean(x, na.rm = TRUE, ...)

# ---- A2. Document, check, install (run inside borgworld project) ----
# In RStudio: Session > Set Working Directory > To Project Directory (borgworld)
# Then run:
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::document()    # regenerates NAMESPACE + Rd
  # devtools::check()     # optional
  devtools::install()     # installs into your user library
} else {
  install.packages("devtools")
}

# ---- A3. Quick reload during development (optional) ----
# Loads package without installing. Run from borgworld project.
if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all()
} else {
  install.packages("pkgload")
}

# ---- A4. Version bump (manual step) ----
# Edit DESCRIPTION and increment Version:. Save file.
# Example: Version: 0.0.0.9001

# =========================================================
# B. Use in a client project
# =========================================================

# ---- B1. Restart R (recommended) ----
# RStudio: Session > Restart R

# ---- B2. Install from source path (one of these) ----
# If borgworld is not installed or you changed it:
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::install(borgworld_path)   # installs current source
}

# Alternative:
# if (requireNamespace("remotes", quietly = TRUE)) {
#   remotes::install_local(borgworld_path)
# }

# renv users:
# if (requireNamespace("renv", quietly = TRUE)) {
#   renv::install(borgworld_path)
#   renv::snapshot()
# }

# ---- B3. Load and use ----
library(borgworld)
# Example usage:
x <- c(1, 2, NA, 4)
bsum(x)    # 7
bmean(x)   # 2.333333

# ---- B4. Notes ----
# - Do not edit NAMESPACE by hand. devtools::document() manages it.
# - Keep borgworld and client projects in separate R sessions if you want
#   to run long tasks in one while working in the other.


# Optional: Additional useful commands

# Check specific aspects
devtools::check_man()  # Check documentation only
devtools::check_examples()  # Run examples only

# If you're using GitHub
devtools::install_github("yourusername/borgworld")  # Install from GitHub
usethis::use_version()  # Bump version number

# Load for interactive testing without installing
devtools::load_all()  # or Ctrl+Shift+L in RStudio

# Check that your examples work
devtools::run_examples()
