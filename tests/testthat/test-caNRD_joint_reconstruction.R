sim_jr <- function(n = 1500, beta = c(0.05, 0, 0), seed = 1) {
  set.seed(seed)
  p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C"))); p <- p / rowSums(p)
  g <- stats::rbinom(n, 2, 0.3)
  Z <- sapply(1:3, function(h) 0.12 + beta[h] * g + stats::rnorm(n, 0, 0.02))
  cv <- stats::rpois(n, 80) + 10
  list(bulk = matrix(stats::rbinom(n, cv, rowSums(p * pmin(pmax(Z, 0), 1))) / cv, 1, n, dimnames = list("site1", rownames(p))),
       g = matrix(g, 1, n, dimnames = list("var1", rownames(p))), p = p, Z = Z,
       theta = matrix(1, 1, 3, dimnames = list("site1", colnames(p))),
       cov = matrix(cv, 1, n, dimnames = list("site1", rownames(p))))
}
slope <- function(z, g) stats::coef(stats::lm(z ~ g))[2]

test_that("caNRD_joint_reconstruction returns documented structure and follows the conditional formula", {
  d <- sim_jr()
  fit <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  rec <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, theta_floor = 0, fit = fit, coverage = d$cov)
  expect_named(rec, c("expected", "reconstructed", "conditional_sd", "residual", "total_sd", "ci_low", "ci_high", "diagnostics"))
  expect_null(rec$total_sd)
  expect_equal(rec$diagnostics$status, "reconstructed")
  expect_equal(dim(rec$reconstructed$A), c(1L, 1500L))
  expect_false(anyNA(rec$reconstructed$B))
  # formula check against a direct computation
  g <- d$g[1, ]; M <- sweep(outer(g, fit$beta), 2, fit$mu, `+`); mb <- rowSums(d$p * M); cv <- d$cov[1, ]
  mc <- pmin(pmax(mb, 0.5 / cv), 1 - 0.5 / cv); v <- mc * (1 - mc) / cv + fit$tau2_0[1]
  A <- sweep(d$p, 2, fit$sigma2, `*`); Zh <- M + sweep(A, 1, (d$bulk[1, ] - mb) / (rowSums(A * d$p) + v), `*`)
  expect_equal(unname(rec$reconstructed$A[1, ]), unname(Zh[, 1]), tolerance = 1e-12)
  expect_equal(unname(rec$expected$C[1, ]), unname(M[, 3]), tolerance = 1e-12)
  expect_true(all(rec$conditional_sd$A >= 0))
})

test_that("known betas remove leakage that a genotype-blind reconstruction shows", {
  d <- sim_jr(seed = 4, beta = c(0.08, 0, 0))
  fit <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  g <- d$g[1, ]
  informed <- fit; informed$beta <- c(0.08, 0, 0)                            # oracle betas, fitted everything else
  # genotype-blind: the same model fitted WITHOUT genotype (mean and variances), beta = 0
  cv <- d$cov[1, ]; fb <- caEditR:::.canrd_eqtl_fit(d$bulk[1, ], d$p, d$p, cv, 0.5 / cv)
  blind <- fit; blind$beta <- 0; blind$mu <- unname(fb$b); blind$sigma2 <- fb$sigma2; blind$tau2_0 <- fb$tau2_0
  ri <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, theta_floor = 0, fit = informed, coverage = d$cov)
  rb <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, theta_floor = 0, fit = blind, coverage = d$cov)
  expect_lt(abs(slope(ri$reconstructed$B[1, ], g)), 0.003)                  # unaffected: no genetic slope
  expect_lt(abs(slope(ri$reconstructed$C[1, ], g)), 0.003)
  expect_gt(slope(ri$reconstructed$A[1, ], g), 0.07)                         # affected: effect kept
  # blind: some unaffected cell type carries a sizeable share of A's effect (which one depends on the fitted variances)
  lb <- max(slope(rb$reconstructed$B[1, ], g), slope(rb$reconstructed$C[1, ], g)) / slope(rb$reconstructed$A[1, ], g)
  expect_gt(lb, 0.1)
  li <- max(abs(slope(ri$reconstructed$B[1, ], g)), abs(slope(ri$reconstructed$C[1, ], g))) / slope(ri$reconstructed$A[1, ], g)
  expect_lt(li, lb / 5)
})

test_that("caNRD_joint_reconstruction input checks and statuses", {
  d <- sim_jr(n = 400)
  fit <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  expect_error(caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, fit = fit), "theta_floor")
  expect_error(caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = data.frame(x = 1)), "caNRD_editQTL")
  two <- rbind(fit, transform(fit, variant_id = "var2"))
  expect_error(caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = two), "one variant per site")
  th <- d$theta; th[, "C"] <- 1e-3
  r <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, th, theta_floor = 1e-3, fit = fit, coverage = d$cov)
  expect_equal(r$diagnostics$status, "gating_mismatch")
  expect_true(all(is.na(r$reconstructed$A)))
})

test_that("bootstrap intervals widen the conditional ones; pooled sigma2 keeps the means", {
  d <- sim_jr(n = 600, seed = 7)
  fit <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  bf <- caNRD_editQTL_bootstrap(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov, n_boot = 8, seed = 1)
  expect_length(bf, 8)
  rb <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = fit, coverage = d$cov, boot_fits = bf)
  expect_equal(rb$diagnostics$n_boot_used, 8L)
  expect_gt(mean(rb$total_sd$A), mean(rb$conditional_sd$A))
  expect_true(all(rb$ci_low$B <= rb$reconstructed$B & rb$ci_high$B >= rb$reconstructed$B))
  rp <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = fit, coverage = d$cov, sigma2 = "pooled")
  expect_equal(rp$expected$A, rb$expected$A)
  expect_error(caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = fit, boot_fits = fit), "boot_fits")
})

test_that("vectorised reconstruction equals the per-site reference implementation", {
  d <- sim_jr(n = 500, seed = 21)
  fit <- caNRD_editQTL(d$bulk, d$g, d$p, d$theta, theta_floor = 0, coverage = d$cov)
  a <- caNRD_joint_reconstruction(d$bulk, d$g, d$p, d$theta, 0, fit = fit, coverage = d$cov, chunk_size = 1)
  b <- caEditR:::.caNRD_joint_reconstruction_reference(d$bulk, d$g, d$p, d$theta, 0, fit = fit, coverage = d$cov)
  for (el in c("expected", "reconstructed", "conditional_sd")) for (h in names(b[[el]])) expect_equal(a[[el]][[h]], b[[el]][[h]], tolerance = 1e-12)
  expect_equal(a$residual, b$residual, tolerance = 1e-12)
  expect_identical(a$diagnostics$status, b$diagnostics$status)
})

test_that("caNRD_editQTL_shrink shrinks null betas more than true effects and keeps the input fit", {
  set.seed(5); S <- 150; n <- 400
  p <- matrix(stats::rgamma(n * 3, 5), n, 3, dimnames = list(paste0("s", 1:n), c("A", "B", "C"))); p <- p / rowSums(p)
  sites <- paste0("site", 1:S); eff <- rep(c(0.06, 0, 0), S); eff[seq(1, 3 * S, 3)][(S / 2 + 1):S] <- 0      # half the sites: effect in A
  g <- matrix(stats::rbinom(S * n, 2, 0.3), S, n, dimnames = list(paste0("v", 1:S), rownames(p)))
  cov <- matrix(stats::rpois(S * n, 80) + 10, S, n, dimnames = list(sites, rownames(p)))
  bulk <- t(sapply(1:S, function(k) { Z <- sapply(1:3, function(h) 0.12 + eff[3 * (k - 1) + h] * g[k, ] + stats::rnorm(n, 0, 0.02))
    stats::rbinom(n, cov[k, ], rowSums(p * pmin(pmax(Z, 0), 1))) / cov[k, ] })); dimnames(bulk) <- list(sites, rownames(p))
  th <- matrix(1, S, 3, dimnames = list(sites, colnames(p)))
  fit <- caNRD_editQTL(bulk, g, p, th, theta_floor = 0, pairs = data.frame(site_id = sites, variant_id = paste0("v", 1:S)), coverage = cov)
  sf <- caNRD_editQTL_shrink(fit, bulk, g, p, th, theta_floor = 0, coverage = cov, min_sites = 50)
  expect_equal(sf$beta_unshrunk, fit$beta)
  expect_true(all(c("pi0", "t2") %in% names(attr(sf, "prior"))))
  trueA <- sf$celltype == "A" & sf$site_id %in% sites[1:(S / 2)]; nullB <- sf$celltype == "B"
  expect_gt(stats::median(abs(sf$beta[trueA] / sf$beta_unshrunk[trueA])), stats::median(abs(sf$beta[nullB] / sf$beta_unshrunk[nullB])))
  expect_lt(mean(abs(sf$beta[nullB])), mean(abs(fit$beta[nullB])))
  rec <- caNRD_joint_reconstruction(bulk, g, p, th, 0, fit = sf, coverage = cov)
  expect_true(all(rec$diagnostics$status == "reconstructed"))
})
