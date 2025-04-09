package_install <- function(packages = kaiserschmRn::dependencies, verbose = 1) {

  stopifnot(is.character(packages), verbose %in% c(0, 1, 2))

  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)

  if (length(to_install) > 0) {

    missing_packages <- to_install
    repo_packages <- rownames(available.packages())
    available_packages <- missing_packages[missing_packages %in% repo_packages]
    unavailable_packages <- setdiff(missing_packages, available_packages)

    if (length(missing_packages) > 0 && length(available_packages) == 0) {
      stop("Packages not available!")
    }

    n_missing <- length(available_packages)
    if (verbose > 0) {
      cat(sprintf("\n\nMissing packages: %s\n", paste(available_packages, collapse = ", ")))
      cat("\nInstalling ...\n")
    }

    for (i in seq_len(n_missing)) {
      if (verbose > 0) {
        cat(sprintf("%0*d of %d: %s\n", nchar(n_missing), i, n_missing, available_packages[i]))
      }

      quiet <- if (verbose == 2) FALSE else TRUE
      install.packages(available_packages[i], quiet = quiet)
    }

    installed <- rownames(installed.packages())
    missing_after <- setdiff(packages, installed)
    if (length(missing_after) == 0) {
      if (verbose > 0) cat("\nDone.\n\n")
    } else {
      warning(sprintf("\n\nPackages not able to install: %s\n\n", paste(missing_after, collapse = ", ")))
    }
  }
  invisible(NULL)
}
