package_install <- function(packages = kaiserschmRn::dependencies, verbose = 1) {

  stopifnot(is.character(packages), verbose %in% c(0, 1, 2))
  github_packages <- grepl("/", packages)
  package_names <- ifelse(github_packages, sub(".*/", "", packages), packages)

  installed <- rownames(installed.packages())
  not_installed_idx <- which(!package_names %in% installed)

  if (length(not_installed_idx) == 0) {
    invisible(NULL)
    return()
  }

  cran_idx <- not_installed_idx[ !grepl("/", packages[not_installed_idx]) ]
  github_idx <- not_installed_idx[ grepl("/", packages[not_installed_idx]) ]

  install_cran_packages <- function(pkg_names) {
    repo_packages <- rownames(available.packages())
    available_pkgs <- pkg_names[pkg_names %in% repo_packages]
    if (length(available_pkgs) == 0) return()

    if (verbose > 0) {
      cat(sprintf("\n\nMissing CRAN packages: %s\n", paste(available_pkgs, collapse = ", ")))
      cat("\nInstalling CRAN packages...\n")
    }
    for (i in seq_along(available_pkgs)) {
      if (verbose > 0) {
        cat(sprintf("%0*d of %d: %s\n", nchar(length(available_pkgs)), i, length(available_pkgs), available_pkgs[i]))
      }
      quiet <- if (verbose == 2) FALSE else TRUE
      install.packages(available_pkgs[i], quiet = quiet)
    }
  }

  install_github_packages <- function(pkgs) {
    if (length(pkgs) == 0) return()

    if (verbose > 0) {
      cat(sprintf("\n\nMissing GitHub packages: %s\n", paste(pkgs, collapse = ", ")))
      cat("\nInstalling GitHub packages...\n")
    }
    for (i in seq_along(pkgs)) {
      if (verbose > 0) {
        cat(sprintf("%0*d of %d: %s\n", nchar(length(pkgs)), i, length(pkgs), pkgs[i]))
      }
      quiet <- if (verbose == 2) FALSE else TRUE
      devtools::install_github(pkgs[i], quiet = quiet)
    }
  }

  if (length(cran_idx) > 0) {
    cran_pkg_names <- package_names[cran_idx]
    install_cran_packages(cran_pkg_names)
  }

  if (length(github_idx) > 0) {
    github_pkgs <- packages[github_idx]  # Originalstring, z.B. "Autor/Paketname"
    install_github_packages(github_pkgs)
  }

  installed <- rownames(installed.packages())
  missing_after <- setdiff(package_names, installed)
  if (length(missing_after) == 0) {
    if (verbose > 0) cat("\nDone.\n\n")
  } else {
    warning(sprintf("\n\nPackages not able to install: %s\n\n", paste(missing_after, collapse = ", ")))
  }
  invisible(NULL)
}
