sim_eq <- function(n = 600, beta = c(0.05, 0, 0), seed = 1) {
  set.seed(seed)
  p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C"))); p <- p / rowSums(p)
  g <- stats::rbinom(n, 2, 0.3)
  Z <- sapply(1:3, function(h) 0.15 + beta[h] * g + stats::rnorm(n, 0, 0.02))
  cv <- stats::rpois(n, 80) + 10
  list(bulk = matrix(stats::rbinom(n, cv, rowSums(p * pmin(pmax(Z, 0), 1))) / cv, 1, n, dimnames = list("site1", rownames(p))),
       g = matrix(g, 1, n, dimnames = list("var1", rownames(p))), p = p,
       theta = matrix(1, 1, 3, dimnames = list("site1", colnames(p))),
       cov = matrix(cv, 1, n, dimnames = list("site1", rownames(p))))
}

test_that("caNRD_editQTL recovers a single-cell-type effect and returns documented columns", {
  d <- sim_eq()
  r <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_equal(nrow(r), 3)
  expect_true(all(c("beta", "se", "ci_low", "ci_high", "p", "p_wald", "p_site", "sigma2", "tau2_0", "vif",
                    "weakly_identifiable", "converged") %in% names(r)))
  expect_true(all(r$status == "tested"))
  expect_true(all(r$converged))
  expect_lt(r$p[r$celltype == "A"], 0.01)
  expect_true(r$ci_low[1] < 0.05 && r$ci_high[1] > 0.05)
  expect_lt(r$p_site[1], 0.01)
  expect_true(all(r$p >= 0 & r$p <= 1))
})

test_that("caNRD_editQTL agrees with celltype_edqtl on the mean model", {
  d <- sim_eq(seed = 3)
  a <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  b <- celltype_edqtl(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_gt(stats::cor(a$beta, b$beta), 0.99)
})

test_that("caNRD_editQTL gating and input checks", {
  d <- sim_eq()
  expect_error(caNRD_editQTL(d$bulk, d$g, d$p, d$theta), "theta_floor")
  th <- d$theta; th[, "C"] <- 1e-3
  r <- caNRD_editQTL(d$bulk, d$g, d$p, th, theta_floor = 1e-3, coverage = d$cov)
  expect_equal(r$status[r$celltype == "C"], "not_identifiable")
  g0 <- d$g; g0[] <- 1
  expect_true(all(caNRD_editQTL(d$bulk, g0, d$p, d$theta, theta_floor = 0)$status == "monomorphic_variant"))
  p2 <- d$p; p2[, "B"] <- p2[, "A"]; p2 <- p2 / rowSums(p2)
  expect_true("aliased" %in% caNRD_editQTL(d$bulk, d$g, p2, d$theta, theta_floor = 0, coverage = d$cov)$status)
})
