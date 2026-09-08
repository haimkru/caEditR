# caEditR

MIT License (see `LICENSE`). Written by Haim Krupkin+Claude.

## Install from GitHub

```bash
git clone https://github.com/haimkru/caEditR.git
```

## Install in R

```r
setwd("~/Downloads")   # the directory that CONTAINS the caEditR folder from the git clone above

for (dep in c("jsonlite", "nnls")) {
  if (!requireNamespace(dep, quietly = TRUE)) {
    install.packages(dep, repos = "https://cloud.r-project.org")
  }
}

install.packages("caEditR", repos = NULL, type = "source")

library(caEditR)
```

## Example

`vignettes/caEditR.Rmd` — a complete worked example. After installing, open it with:

```r
vignette("caEditR")
```
