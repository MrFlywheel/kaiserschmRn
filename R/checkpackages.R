dependencies <- c('lidR', 'terra', 'feather', 'piledge/pizzR', "parallel", 'doParallel', 'foreach')

to_install <- !dependencies %in% rownames(installed.packages())

if (any(to_install)) {
  missing_dependencies <- as.character(dependencies[to_install])
  cat("\n\nMissing dependencies:", paste(missing_dependencies, collapse = ", "), "\n")
  cat("Use 'kaiserschmaRn::package_install(kaiserschRn::dependencies)' to add missing packages.\n\n")
} else {
  cat("All dependencies are installed.\n")
}
