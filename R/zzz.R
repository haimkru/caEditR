#' @keywords internal
"_PACKAGE"

.caEditR_env <- new.env(parent = emptyenv())

#' A CRAN mirror known to actually work, used ONLY for auto-installs.
#'
#' `getOption("repos")` cannot be trusted -- it's frequently left pointing
#' at a broken/blocked mirror by a user's own `.Rprofile` (seen directly in
#' this project's own dev environment: a hardcoded plain-http
#' `cran.us.r-project.org` that 403s). `.ensure_installed()` temporarily
#' overrides `options(repos=...)` to this mirror for the duration of an
#' automatic install, then restores whatever was there before -- so
#' auto-install works regardless of what the calling session's own
#' `.Rprofile`/environment has configured.
#' @keywords internal
.CRAN_MIRROR <- "https://cloud.r-project.org"

#' Run `expr` with `options(repos=...)` temporarily forced to a working
#' CRAN mirror (`.CRAN_MIRROR`), restoring the previous value afterward --
#' used by `.ensure_installed()` so `install.packages()`/`BiocManager::install()`/
#' `remotes::install_github()` (all of which read `options("repos")`) work
#' even when the ambient repos option is broken.
#'
#' IMPORTANT: when `BiocManager` is available, this overrides ONLY the
#' `CRAN` entry within `BiocManager::repositories()`, not the whole repos
#' vector -- replacing it entirely (as an earlier version of this function
#' did) breaks `BiocManager::install()`'s ability to resolve Bioconductor
#' -only transitive dependencies (e.g. installing `GenomicFeatures` needs
#' `Rhtslib`/`Rsamtools`/`rtracklayer`, which only exist in Bioconductor's
#' own repos, not CRAN's) -- confirmed directly: replacing the whole repos
#' vector made exactly that install fail with "package 'Rhtslib' is not
#' available for this version of R".
#' @keywords internal
.with_reliable_cran <- function(expr) {
  old <- options()
  on.exit(options(old), add = TRUE)
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    repos <- BiocManager::repositories()
    repos["CRAN"] <- .CRAN_MIRROR
    options(repos = repos)
  } else {
    options(repos = c(CRAN = .CRAN_MIRROR))
  }
  force(expr)
}

#' Locate a Python interpreter with numpy/scipy, without using reticulate.
#'
#' Tries, in order: (1) `options(caEditR.python=...)`, (2) the
#' `CAEDITR_PYTHON` environment variable, (3) `~/envs/catca-edit/bin/python3`
#' if it exists (this project's own dev environment), (4) whatever
#' `python3` is first on PATH.
#'
#' Deliberately does NOT use `reticulate::conda_list()`/`reticulate::use_python()`
#' -- this package no longer embeds Python inside the R process at all (see
#' `.run_python_op()`'s own docstring for why), so there is nothing to
#' "configure" beyond finding a path to hand to `system2()`.
#' @keywords internal
.find_python <- function() {
  candidate <- getOption("caEditR.python", Sys.getenv("CAEDITR_PYTHON", NA))
  if (is.na(candidate) || !nzchar(candidate)) {
    guess <- path.expand("~/envs/catca-edit/bin/python3")
    if (file.exists(guess)) candidate <- guess
  }
  if (is.na(candidate) || !nzchar(candidate)) {
    candidate <- Sys.which("python3")
    if (!nzchar(candidate)) candidate <- Sys.which("python")
  }
  if (is.na(candidate) || !nzchar(candidate)) {
    stop("No Python interpreter found. Set options(caEditR.python=\"/path/to/python3\") ",
         "or Sys.setenv(CAEDITR_PYTHON=\"/path/to/python3\") to a Python with numpy and scipy installed.",
         call. = FALSE)
  }
  candidate
}

#' The `lib/` directory of the Python's own conda/venv environment, if any.
#' @keywords internal
.python_env_lib_dir <- function(python_path) {
  env_root <- dirname(dirname(normalizePath(python_path, mustWork = FALSE)))  # <env_root>/bin/python3 -> <env_root>
  lib_dir <- file.path(env_root, "lib")
  if (dir.exists(lib_dir)) lib_dir else NA_character_
}

#' Ensure an R package is installed, installing it automatically if not.
#'
#' Used by `TCA_Like()` (for the real CRAN `TCA` package) and
#' `build_music_reference()`/`estimate_proportions_music()` (for `MuSiC`
#' and its Bioconductor dependencies) so that a first-time call "just
#' works" without the user having to separately run `install.packages()`
#' themselves first -- per direct user request. Always prints a clear
#' `message()` before installing anything (never silent), and raises an
#' informative error (not a cryptic one) if the automatic install itself
#' fails (e.g. no network access) or still doesn't produce a loadable
#' package.
#'
#' @param pkg package name to check/install.
#' @param install_fn a zero-argument function that installs `pkg` if called.
#' @keywords internal
.ensure_installed <- function(pkg, install_fn) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  message(sprintf(
    "caEditR: the '%s' package is required but not installed. Installing it automatically now (per this package's design -- see ?%s)...",
    pkg, pkg
  ))
  result <- tryCatch({ .with_reliable_cran(install_fn()); TRUE }, error = function(e) e)
  # A failed compiled-code install (e.g. a package with C/C++ source, like
  # MuSiC's MCMCpack dependency, or GenomicFeatures' Rhtslib) is one of the
  # few auto-install failure modes that ISN'T fixable by anything caEditR's
  # own code can do -- it means the R installation's own compiler toolchain
  # is broken, which affects every package trying to compile code there,
  # not just this one. Pointing at this project's own pinned conda
  # environment (verified, this whole package's own dev/test environment)
  # is the most actionable thing to say in that situation.
  toolchain_hint <- paste(
    "\nIf this looks like a compiler/toolchain failure (e.g. mentions",
    "'configure: error', 'compilation failed', or similar above) rather",
    "than a missing package or network issue: that means this R",
    "installation's own C/C++ compiler setup is broken, which no R",
    "package's code can fix. Try this project's own pinned, tested conda",
    "environment instead (see README.md's 'Install' section):",
    "conda env create -f environment.yml && conda activate catca-edit."
  )
  if (inherits(result, "error")) {
    stop(sprintf(
      "caEditR: automatic installation of '%s' failed (%s). Install it manually and retry.%s",
      pkg, conditionMessage(result), toolchain_hint
    ), call. = FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf(
      "caEditR: '%s' still not available after attempting automatic installation. Install it manually and retry.%s",
      pkg, toolchain_hint
    ), call. = FALSE)
  }
  message(sprintf("caEditR: '%s' installed successfully.", pkg))
  invisible(TRUE)
}

#' Run one operation of the vendored Python math in a FRESH subprocess.
#'
#' THE CORE DESIGN CHOICE THIS FUNCTION EXISTS TO IMPLEMENT: earlier
#' versions of this package called the vendored core.py/no_reference.py/
#' simulate.py IN-PROCESS via `reticulate` (Python embedded directly
#' inside the running R session). That is simpler when it works, but
#' fails unpredictably in a long-running R session (e.g. RStudio Server
#' left open for a while): R's OWN process may have already loaded an
#' incompatible system `libstdc++.so.6` for entirely unrelated reasons
#' (its own C++ dependencies, or another already-loaded R package) before
#' this package's code ever runs. Once a shared library with a given name
#' is mapped into a process, the dynamic linker reuses that SAME mapping
#' for everything else requesting it (e.g. Python's compiled scipy
#' extension) -- no `Sys.setenv(LD_LIBRARY_PATH=...)` from anywhere in
#' that already-running process can undo this. This is a well-documented,
#' repeatedly-reported reticulate failure mode (rstudio/reticulate issues
#' #1467, #841, #311, #428), not something fixable purely from R code
#' running inside the affected process -- confirmed directly during this
#' package's own development (the in-process version worked in a fresh
#' `Rscript` call but failed in a long-running interactive session with
#' exactly this class of error).
#'
#' THE FIX: never embed Python inside R's process. Every exported
#' function that needs this vendored math spawns a BRAND NEW, short-lived
#' Python subprocess (`inst/python/cli_driver.py`) via `system2()`, with
#' its OWN environment (including `LD_LIBRARY_PATH`) set explicitly at
#' THAT subprocess's own exec() time -- a genuinely fresh OS process,
#' completely unaffected by whatever R's own process has already loaded.
#' This fixes the root cause rather than working around a symptom, and
#' means `library(caEditR)` plus a direct function call works the same
#' way regardless of how long the R session has been running or what
#' else has been loaded into it.
#'
#' @param op operation name (matches a key in `cli_driver.py`'s `_OPS` dict).
#' @param args named list of JSON-serializable arguments for that operation.
#' @return the parsed JSON response (as an R list/vector/matrix via
#'   `jsonlite::fromJSON`'s automatic simplification).
#' @keywords internal
.run_python_op <- function(op, args) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required. Install with install.packages('jsonlite').", call. = FALSE)
  }
  python_path <- .find_python()
  lib_dir <- .python_env_lib_dir(python_path)
  driver <- system.file("python", "cli_driver.py", package = "caEditR")
  if (!nzchar(driver)) stop("Could not locate inst/python/cli_driver.py in the installed caEditR package.", call. = FALSE)

  input_file <- tempfile(fileext = ".json")
  output_file <- tempfile(fileext = ".json")
  on.exit(unlink(c(input_file, output_file)), add = TRUE)
  jsonlite::write_json(args, input_file, auto_unbox = TRUE, digits = NA, na = "null", null = "null")

  env <- character(0)
  if (!is.na(lib_dir)) {
    current_ld <- Sys.getenv("LD_LIBRARY_PATH")
    env <- paste0("LD_LIBRARY_PATH=", lib_dir, if (nzchar(current_ld)) paste0(":", current_ld) else "")
  }

  out <- suppressWarnings(system2(
    python_path, args = shQuote(c(driver, op, input_file, output_file)),
    env = env, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0 && !file.exists(output_file)) {
    stop("caEditR's Python subprocess (", python_path, ") exited with status ", status,
         " and produced no output.\nCaptured output:\n", paste(out, collapse = "\n"), call. = FALSE)
  }
  if (!file.exists(output_file)) {
    stop("caEditR's Python subprocess produced no output file.\nCaptured output:\n",
         paste(out, collapse = "\n"), call. = FALSE)
  }
  result <- jsonlite::fromJSON(output_file, simplifyVector = TRUE)
  if (!is.null(result$error)) {
    stop("caEditR's Python subprocess raised an error:\n", result$error, call. = FALSE)
  }
  result
}
