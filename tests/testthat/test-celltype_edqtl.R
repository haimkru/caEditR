sim_edqtl_data <- function(n = 600, beta_target = 0.05, seed = 1, zero_donors = 0) {
  set.seed(seed)
  cts <- c("Mono", "Neut", "NK")
  p <- matrix(stats::rgamma(n * 3, shape = c(3, 5, 3)), n, 3, byrow = TRUE, dimnames = list(paste0("s", seq_len(n)), cts))
  p <- p / rowSums(p)
  if (zero_donors > 0) { p[seq_len(zero_donors), c("Mono", "Neut")] <- 0; p[seq_len(zero_donors), "NK"] <- 1 }
  theta <- matrix(c(3, 2, 2), 1, 3, dimnames = list("site1", cts))
  g <- stats::rbinom(n, 2, 0.35)
  e <- cbind(Mono = 0.10 + beta_target * g, Neut = rep(0.08, n), NK = rep(0.12, n))
  e <- pmin(pmax(e + matrix(stats::rnorm(n * 3, 0, 0.02), n, 3), 0), 1)
  w <- sweep(p, 2, theta[1, ], `*`); phi <- w / rowSums(w)
  cov <- stats::rpois(n, 60) + 10
  k <- stats::rbinom(n, cov, rowSums(phi * e))
  list(bulk = matrix(k / cov, 1, n, dimnames = list("site1", rownames(p))),
       cov = matrix(cov, 1, n, dimnames = list("site1", rownames(p))),
       geno = matrix(g, 1, n, dimnames = list("var1", rownames(p))), p = p, theta = theta)
}

test_that("celltype_edqtl finds a one-cell-type effect in that cell type only, with the right size", {
  d <- sim_edqtl_data(beta_target = 0.05)
  r <- celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_equal(nrow(r), 3)
  expect_true(all(r$status == "tested"))
  mono <- r[r$celltype == "Mono", ]
  expect_lt(mono$p, 1e-6)
  expect_equal(mono$beta, 0.05, tolerance = 0.3)  # relative tolerance
  expect_true(all(r$p[r$celltype != "Mono"] > 1e-3))
})

test_that("celltype_edqtl keeps false positives controlled when there is no effect", {
  ps <- unlist(lapply(1:40, function(s) {
    d <- sim_edqtl_data(beta_target = 0, seed = 100 + s)
    celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 0, coverage = d$cov)$p
  }))
  expect_lt(mean(ps < 0.05), 0.12)
})

test_that("celltype_edqtl gates floored / low-phi cell types like caNRD_edit and requires theta_floor", {
  d <- sim_edqtl_data()
  d$theta[1, "NK"] <- 1e-3
  r <- celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 1e-3, coverage = d$cov)
  expect_equal(r$status[r$celltype == "NK"], "not_identifiable")
  expect_true(all(is.na(r[r$celltype == "NK", c("beta", "p")])))
  expect_equal(unique(r$n_celltypes), 2)
  expect_error(celltype_edqtl(d$bulk, d$geno, d$p, d$theta), "theta_floor is required")
})

test_that("celltype_edqtl drops donors with zero proportion for all identifiable cell types", {
  d <- sim_edqtl_data(zero_donors = 15)
  d$theta[1, "NK"] <- 1e-3
  r <- celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 1e-3, coverage = d$cov)
  expect_equal(unique(r$n_samples), 600 - 15)
  expect_true(all(is.finite(r$p[r$status == "tested"])))
})

test_that("celltype_edqtl labels degenerate sites instead of returning NaN p-values", {
  d <- sim_edqtl_data()
  b <- d$bulk; b[1, ] <- 0
  r <- celltype_edqtl(b, d$geno, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_true(all(r$status == "no_variation_in_bulk"))
  expect_true(all(is.na(r$p)))
})

test_that("celltype_edqtl flags perfectly collinear cell types as aliased", {
  d <- sim_edqtl_data()
  d$p[, "Neut"] <- d$p[, "Mono"]; d$p <- d$p / rowSums(d$p); d$theta[1, "Neut"] <- d$theta[1, "Mono"]
  r <- celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_true(any(r$status == "aliased"))
  expect_true(all(is.na(r$p[r$status == "aliased"])))
})

test_that("celltype_edqtl drops samples with missing covariates and reports convergence", {
  d <- sim_edqtl_data()
  covs <- matrix(stats::rnorm(600), 600, 1, dimnames = list(colnames(d$bulk), "pc1")); covs[1:7, 1] <- NA
  r <- celltype_edqtl(d$bulk, d$geno, d$p, d$theta, theta_floor = 0, coverage = d$cov, covariates = covs)
  expect_equal(unique(r$n_samples), 593)
  expect_true(all(r$converged[r$status == "tested"]))
})
