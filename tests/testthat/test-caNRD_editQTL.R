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

test_that("regressions from the adversarial review (09.28.2026)", {
  d <- sim_eq(seed = 5)
  a <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  # coverage = NULL: TCA's variance model, SEs comparable to the coverage-aware fit (previously ~8x too large)
  b <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0)
  expect_lt(max(abs(log(b$se / a$se))), log(1.5))
  # unnamed init in cell-type order works; named init is equivalent
  u <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov, init = list(sigma2 = c(4e-4, 4e-4, 4e-4)))
  expect_true(all(u$status == "tested"))
  expect_equal(u$beta, a$beta, tolerance = 1e-3)
  expect_error(caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, init = list(sigma2 = c(1, 1))), "init")
  # coverage without dimnames: informative error
  expect_error(caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = unname(d$cov)), "coverage must have")
  # near-collinear cell types are aliased, not tested with huge betas
  p3 <- d$p; p3[, "B"] <- 0.999 * p3[, "A"] + 0.001 * p3[, "B"]; p3 <- p3 / rowSums(p3)
  r3 <- caNRD_editQTL(d$bulk, d$g, p3, d$theta, theta_floor = 0, coverage = d$cov)
  expect_true(all(r3$status[r3$celltype %in% c("A", "B")] == "aliased"))
  # LRT statistics are non-negative by construction (nesting guard); p == 1 only if beta is exactly 0
  expect_true(all(a$p < 1))
})

test_that("fast engine matches the reference engine", {
  d <- sim_eq(seed = 11, beta = c(0.04, 0, 0.02))
  f <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  r <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov, engine = "reference")
  expect_identical(f$status, r$status)
  expect_equal(f$beta, r$beta, tolerance = 1e-3)
  expect_equal(f$se, r$se, tolerance = 1e-3)
  expect_gt(f$loglik[1], r$loglik[1] - 1e-3)
  expect_error(caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, engine = "reference", chunk_size = 5), "fast")
})

test_that("scan engine agrees with the exact engine", {
  set.seed(31); n <- 500; S <- 6
  p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C"))); p <- p / rowSums(p)
  sites <- paste0("site", 1:S)
  G <- matrix(stats::rbinom(4 * S * n, 2, 0.3), 4 * S, n, dimnames = list(paste0("v", 1:(4 * S)), rownames(p)))
  cov <- matrix(stats::rpois(S * n, 70) + 10, S, n, dimnames = list(sites, rownames(p)))
  bulk <- t(sapply(1:S, function(k) { g <- G[4 * k - 3, ]; Z <- sapply(1:3, function(h) 0.12 + (h == 1) * 0.04 * g + stats::rnorm(n, 0, 0.02))
    stats::rbinom(n, cov[k, ], rowSums(p * pmin(pmax(Z, 0), 1))) / cov[k, ] })); dimnames(bulk) <- list(sites, rownames(p))
  th <- matrix(1, S, 3, dimnames = list(sites, colnames(p)))
  pairs <- data.frame(site_id = rep(sites, each = 4), variant_id = paste0("v", 1:(4 * S)))
  ex <- caNRD_editQTL(bulk, G, p, th, 0, pairs = pairs, coverage = cov)
  sc <- caNRD_editQTL(bulk, G, p, th, 0, pairs = pairs, coverage = cov, engine = "scan", refine = NULL)
  sr <- caNRD_editQTL(bulk, G, p, th, 0, pairs = pairs, coverage = cov, engine = "scan", refine = 1e-3)
  sl <- caNRD_editQTL(bulk, G, p, th, 0, pairs = pairs, coverage = cov, engine = "scan")
  expect_equal(sum(sl$refined) / 3, length(unique(sl$site_id[is.finite(sl$p_site)])))
  expect_identical(sc$status, ex$status)
  expect_gt(stats::cor(sc$beta, ex$beta, use = "complete.obs"), 0.99)
  expect_true(all(sr$refined[sr$p < 1e-3 & !is.na(sr$p)]))
  rr <- sr$refined; expect_equal(sr$beta[rr], ex$beta[rr], tolerance = 1e-8)
  expect_false(is.null(attr(caNRD_editQTL(bulk, G, p, th, 0, pairs = pairs, coverage = cov, engine = "scan", vcov = "beta"), "coef_cov")))
})

test_that("scan engine rejects non-finite / out-of-range genotypes like the exact engine", {
  set.seed(2); n <- 200
  p <- matrix(stats::rgamma(n * 2, 5), n, 2, dimnames = list(paste0("s", 1:n), c("A", "B"))); p <- p / rowSums(p)
  g <- matrix(stats::rbinom(2 * n, 2, 0.3), 2, n, dimnames = list(c("v1", "v2"), rownames(p))); g[1, 3] <- Inf
  y <- matrix(stats::runif(n, 0.05, 0.2), 1, n, dimnames = list("site1", rownames(p))); th <- matrix(1, 1, 2, dimnames = list("site1", c("A", "B")))
  expect_error(caNRD_editQTL(y, g, p, th, 0, engine = "scan"), "dosages")
})
