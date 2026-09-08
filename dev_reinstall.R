# Reinstall caEditR from THIS source tree into your current R session's own
# library, then load it fresh -- run this at the start of every session
# while developing against the source directly (instead of a released/
# tarball version), so `library(caEditR)` always reflects the latest
# source, no matter which R installation or library path your session
# happens to be using.
#
# Deliberately explicit about WHICH library gets the install and WHICH
# library gets loaded, and VERIFIES the reload actually picked up the new
# source (not a stale copy) before declaring success -- confirmed directly
# to matter: a plain `install.packages(pkg_dir, repos=NULL, type="source")`
# with no `lib=` argument can silently install into a DIFFERENT library
# than the one `library(caEditR)` subsequently loads from, when something
# earlier in the session (e.g. a user's own ~/.Rprofile prepending an
# unrelated library path) has reordered .libPaths() -- the classic
# "reinstalled it but nothing changed" symptom. This script pins BOTH
# steps to the exact same, explicit library, then fails loudly (not
# silently) if verification still doesn't find the fresh source.
#
# Usage (works on ANY machine, not just this cluster -- see the portable
# self-location logic below):
#   source("/path/to/wherever/you/copied/caEditR/dev_reinstall.R")
# or, from inside the caEditR/ directory itself:
#   source("dev_reinstall.R")

pkg_dir <- (function() {
  # Self-locate THIS script's own file path when source()d -- the most
  # portable option, since it doesn't depend on the working directory or
  # this specific cluster's own paths at all. Standard R trick: when
  # source()d, the enclosing call frame carries the file path in $ofile.
  sourced_from <- tryCatch({
    frames <- sys.frames()
    ofiles <- vapply(frames, function(fr) {
      of <- tryCatch(fr$ofile, error = function(e) NULL)
      if (is.null(of)) NA_character_ else of
    }, character(1))
    ofiles <- ofiles[!is.na(ofiles)]
    if (length(ofiles) > 0) dirname(ofiles[length(ofiles)]) else NA_character_
  }, error = function(e) NA_character_)

  candidates <- c(
    sourced_from,
    ".",  # if you cd'd into caEditR/ first, per the Usage note above
    "/oak/stanford/groups/smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/caEditR",
    "/labs/smontgom/grps_smontgom/hkrupkin/RNA_Editing_deconvolution/catca-edit-hpc/caEditR"
  )
  for (d in candidates) if (!is.na(d) && file.exists(file.path(d, "DESCRIPTION"))) return(normalizePath(d))
  stop("Could not locate the caEditR source directory (looked for a DESCRIPTION file). ",
       "If you copied this package to a new machine, either cd into the ",
       "caEditR/ folder first, or source() this file using its full path.", call. = FALSE)
})()

if ("caEditR" %in% loadedNamespaces()) {
  try(unloadNamespace("caEditR"), silent = TRUE)
}

# Pin the target library EXPLICITLY: reuse wherever caEditR is already
# installed (so we overwrite the exact copy that would otherwise get
# loaded), or fall back to the first library in .libPaths() if it isn't
# installed anywhere yet.
target_lib <- tryCatch(dirname(find.package("caEditR", quiet = TRUE)), error = function(e) NULL)
if (is.null(target_lib) || !length(target_lib)) target_lib <- .libPaths()[1]
if (!dir.exists(target_lib)) dir.create(target_lib, recursive = TRUE)

# Make sure THIS library is searched first when we load below, regardless
# of what a ~/.Rprofile may have prepended ahead of it earlier in the session.
.libPaths(c(target_lib, .libPaths()))

try(remove.packages("caEditR", lib = target_lib), silent = TRUE)
# Pre-install caEditR's own Imports (jsonlite) from a real CRAN mirror
# FIRST, before installing caEditR itself -- a real failure, seen directly
# on a fresh macOS R install with no packages yet ("dependency 'jsonlite'
# is not available"). Passing a real `repos=` straight to the caEditR
# install call itself does NOT fix this -- confirmed directly that once
# `repos` is non-NULL, install.packages() treats `pkgs` as a package NAME
# to look up in that repo instead of a real local directory ("package
# '<path>' is not available for this version of R"). So `repos = NULL` is
# required for the local caEditR install; the fix is ensuring its
# dependencies are already satisfied before that call runs.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
install.packages(pkg_dir, repos = NULL, type = "source", lib = target_lib)

library(caEditR, lib.loc = target_lib)

# Verify the reload actually picked up THIS source tree, not some other
# stale copy still sitting earlier in .libPaths() -- check for a file that
# only exists in the current source (bump this check whenever you add a
# new marker file/function you know is recent).
loaded_from <- find.package("caEditR")
marker <- file.path(loaded_from, "extdata", "gse64655_bulk_editing_ratios.csv")
if (!file.exists(marker)) {
  stop(sprintf(paste(
    "caEditR reinstalled into '%s' and loaded from '%s', but a file",
    "known to exist in the current source (%s) is still missing there.",
    "This almost always means another, older caEditR install is earlier",
    "in .libPaths() and is shadowing the one just installed. Run",
    "'.libPaths()' and 'find.package(\"caEditR\")' to see which libraries",
    "exist and which one is actually being loaded, then either remove the",
    "stale copy or reorder .libPaths() so '%s' comes first."
  ), target_lib, loaded_from, basename(marker), target_lib), call. = FALSE)
}

message("caEditR reinstalled from source and loaded fresh: ", pkg_dir,
        " (library: ", target_lib, ", verified up to date)")

# Install TCA alongside caEditR itself, right here (AFTER install+load has
# fully finished), rather than waiting for the first TCA_Like() call deep
# in a script. Deliberately NOT done via .onAttach()/.onLoad() inside the
# package itself: R's own install.packages()/R CMD INSTALL machinery calls
# library() internally, MORE THAN ONCE, purely to verify the package can be
# loaded ("testing if installed package can be loaded from temporary/final
# location") -- an .onAttach() hook fires during THOSE internal test-loads
# too, which was verified directly to cause a real failure (a staged-install
# path conflict, since the network install landed mid-install-process,
# inside caEditR's own staging directory). Doing it here instead -- after
# `library(caEditR)` above has already fully completed -- avoids that
# entirely while still achieving "install alongside caEditR" for real use.
caEditR:::.ensure_installed("TCA", function() utils::install.packages("TCA"))
