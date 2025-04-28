dependencies <- c('moments', 'lidR', 'terra', 'feather', 'piledge/pizzR', "parallel", 'doParallel', 'foreach', 'geometry')

installed_pkgs <- installed.packages()[, "Package"]
missing_deps  <- setdiff(dependencies, installed_pkgs)

if (length(missing_deps) > 0) {
  col_start <- "\033[31m"
  col_end   <- "\033[39m"
  bold_on  <- "\033[1m"
  bold_off <- "\033[22m"

  cat(sprintf("%s%s%s%s%s", col_start, "\nMissing dependencies: ", paste(missing_deps, collapse = ", "), "\n", col_end))
  cat(sprintf("%s%s%s%s%s%s%s", col_start, "Use '", bold_on, "pizzR::package.install(pizzR::dependencies)", bold_off, "' to add missing packages.\n\n", col_end))
} else {
  cat("All dependencies are installed.\n")
}
