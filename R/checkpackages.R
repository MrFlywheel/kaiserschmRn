dependencies <- c('moments', 'lidR', 'terra', 'feather', 'piledge/pizzR', "parallel", 'doParallel', 'foreach')

installed_pkgs <- installed.packages()[, "Package"]
missing_deps  <- setdiff(dependencies, installed_pkgs)

if (length(missing_deps) > 0) {
  cat("\nMissing dependencies:", paste(missing_deps, collapse = ", "), "\n")
  cat("Use 'kaiserschmRn::package.install(pizzR::dependencies)' to add missing packages.\n\n")
} else {
  cat("All dependencies are installed.\n")
}
