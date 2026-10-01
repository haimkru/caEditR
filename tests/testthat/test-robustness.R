# Regressions found by the pre-release adversarial review (simulated data only).
co <- simulate_edqtl_cohort(200, n_sites_per_type = 2, n_variants = 2, seed = 11)
Y <- co$bulk_editing; CV <- co$coverage; P <- co$proportions; TH <- co$theta

test_that("caNRD_edit() handles NA editing, NA coverage and an all-NA site with both estimators", {
  Y1 <- Y; Y1[1, 1] <- NA; Y1[2, ] <- NA; CV1 <- CV; CV1[3, 4] <- NA
  for (est in c("ml", "moment")) {
    r <- suppressMessages(caNRD_edit(Y1, CV1, P, TH, theta_floor = 0, estimator = est))
    expect_equal(r$diagnostics$status[r$diagnostics$site_id == rownames(Y)[2]], "skipped: fewer usable samples than cell types")
    expect_true(all(is.na(r$deconvolved[[1]][1, 1])))
    expect_true(any(vapply(r$deconvolved, function(m) any(is.finite(m[4, ])), logical(1))))
  }
})

test_that("caNRD_edit() and caRD_edit()-style inputs align coverage columns by sample id", {
  sh <- rev(seq_len(ncol(CV)))
  for (est in c("ml", "moment")) {
    a <- suppressMessages(caNRD_edit(Y, CV, P, TH, theta_floor = 0, estimator = est))
    b <- suppressMessages(caNRD_edit(Y, CV[, sh], P, TH, theta_floor = 0, estimator = est))
    expect_equal(a$deconvolved, b$deconvolved)
    expect_equal(a$diagnostics, b$diagnostics)
  }
  expect_error(caNRD_edit(Y, unname(CV), P, TH, theta_floor = 0), "row names")
})

test_that("an all-NA site gives not_identifiable instead of aborting edQTL runs", {
  Ya <- Y; Ya[1, ] <- NA; s1 <- rownames(Y)[1]
  for (eng in c("scan", "fast")) {
    r <- caNRD_editQTL(Ya, co$genotypes, P, TH, theta_floor = 0, pairs = co$pairs, coverage = CV, engine = eng)
    expect_true(all(r$status[r$site_id == s1] == "not_identifiable"))
  }
  lead <- co$pairs[!duplicated(co$pairs$site_id), ]
  fit <- caNRD_editQTL(Y, co$genotypes, P, TH, theta_floor = 0, pairs = lead, coverage = CV)
  jr <- caNRD_joint_reconstruction(Ya, co$genotypes, P, TH, theta_floor = 0, fit = fit, coverage = CV)
  expect_equal(jr$diagnostics$status[jr$diagnostics$site_id == s1], "not_identifiable")
})

test_that("scan impute_genotypes = 'mean' matches 'none' closely on sparse missingness", {
  G <- co$genotypes; G[cbind(seq_len(nrow(G)), rep(1:3, length.out = nrow(G)))] <- NA
  a <- caNRD_editQTL(Y, G, P, TH, theta_floor = 0, pairs = co$pairs, coverage = CV, engine = "scan", refine = NULL)
  b <- caNRD_editQTL(Y, G, P, TH, theta_floor = 0, pairs = co$pairs, coverage = CV, engine = "scan", refine = NULL, impute_genotypes = "mean")
  expect_gt(cor(a$beta, b$beta, use = "complete.obs"), 0.99)
})

test_that("simulators and bootstrap leave the user's random-number state unchanged", {
  set.seed(99); before <- .Random.seed
  simulate_edqtl_cohort(50, n_sites_per_type = 1, n_variants = 1, seed = 5)
  expect_identical(.Random.seed, before)
  a <- simulate_edqtl_cohort(50, n_sites_per_type = 1, n_variants = 1, seed = 5)
  b <- simulate_edqtl_cohort(50, n_sites_per_type = 1, n_variants = 1, seed = 5)
  expect_identical(a$bulk_editing, b$bulk_editing)
})

test_that("helpers validate their inputs", {
  expect_error(binomial_tau2(c(.1, .2, .3), c(10, 20)), "lengths")
  expect_true(is.na(binomial_tau2(0.2, NA)))
  expect_error(compute_effective_weights(c(.5, .5), c(2, -1)), "non-negative")
})
