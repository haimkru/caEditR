#' Simulate bulk host-gene expression from proportions and per-cell-type expression weights
#'
#' Bulk expression of each site's host gene in each sample: the sample's cell-type proportions times the per-cell-type
#' expression weights theta, with multiplicative log-normal noise. This is the input [estimate_theta_nnls()] expects.
#'
#' @param proportions samples x cell types matrix (rows sum to 1), with row and column names.
#' @param theta sites x cell types matrix of relative expression weights (0 = not expressed), with row names.
#' @param noise_sdlog sd of the log-normal noise (default 0.15).
#' @param seed random seed (default 0).
#' @return sites x samples matrix of bulk expression.
#' @examples
#' p <- simulate_proportions(10, c(0.5, 0.3, 0.2))
#' dimnames(p) <- list(paste0("s", 1:10), c("A", "B", "C"))
#' th <- matrix(c(1, 2, 0, 1, 1, 1), 2, 3, byrow = TRUE, dimnames = list(c("1:100:+", "2:200:-"), colnames(p)))
#' simulate_bulk_expression(p, th)
#' @export
simulate_bulk_expression <- function(proportions, theta, noise_sdlog = 0.15, seed = 0) {
  proportions <- as.matrix(proportions); theta <- as.matrix(theta)
  if (!all(colnames(proportions) %in% colnames(theta))) stop("theta must have a column for every cell type in proportions", call. = FALSE)
  set.seed(seed)
  e <- tcrossprod(theta[, colnames(proportions), drop = FALSE], proportions)
  e * matrix(exp(stats::rnorm(length(e), 0, noise_sdlog)), nrow(e))
}

#' Simulate a cohort for cell-type edQTL analysis (genotypes, bulk editing, coverage, truth)
#'
#' Generates a cohort with known cell-type genetic effects on RNA editing, in the formats every caEditR function
#' uses: sites x donors `bulk_editing` and `coverage`, variants x donors `genotypes`, donors x cell types
#' `proportions`, sites x cell types `theta`, a `pairs` table (site, variant) and `bulk_expression` of the host genes.
#' Site ids follow [format_site_id()].
#'
#' Design (each site has its own variants; one of them is causal):
#' \itemize{
#'   \item site types: an effect in exactly one cell type (`only_<celltype>`, one class per cell type), the same effect
#'     in every cell type that expresses the host gene (`shared`), or no effect (`null`);
#'   \item per site, `n_variants` cis variants: variant 1 is causal (MAF `maf`), variants 2-3 are in LD with it
#'     (each allele copied from the causal haplotype with probability 0.9 and 0.6), the rest are independent
#'     (MAF uniform on 0.05-0.5);
#'   \item latent editing of donor i in cell type h: \eqn{Z_{ih} = \mu_h + \beta_h G_i + \epsilon_{ih}}, clipped to
#'     \[0, 1\]; baselines \eqn{\mu_h} = U(`baseline`) x U(0.6, 1.4) per site and cell type; sd = `sd` x U(0.5, 1.5);
#'   \item host-gene expression weights theta ~ lognormal(0, 0.5); with probability `theta_zero_p` a cell type does
#'     not express the host gene at a site (theta = 0; at least 2 cell types express it; a cell type with an effect
#'     always does). RNA shares \eqn{\phi = p \theta / \sum p \theta};
#'   \item reads: `coverage = "realistic"`: Poisson(site depth x donor sequencing depth x the donor's relative host-gene
#'     expression), site depth lognormal (median `mean_coverage`, log-sd 1, clipped 5-1000), donor depth
#'     lognormal(0, 0.6), at least 2 reads; `"poisson"`: Poisson(`mean_coverage`) + 10. Edited reads
#'     ~ Binomial(coverage, \eqn{\sum_h \phi_{ih} Z_{ih}}).
#' }
#'
#' @param n_donors number of donors (default 1000).
#' @param celltypes cell-type names (default four blood cell types).
#' @param composition mean cell-type proportions (Dirichlet mean), same length as `celltypes`.
#' @param concentration Dirichlet concentration (default 20).
#' @param n_sites_per_type number of sites of each type (default 10): `length(celltypes)` single-cell-type classes,
#'   shared and null.
#' @param n_variants cis variants per site (default 5, at least 1).
#' @param maf allele frequency of the causal variant (default 0.25).
#' @param beta per-allele effect on the editing ratio (default 0.05).
#' @param baseline range of the site's baseline editing (default 0.05-0.30).
#' @param sd biological sd of latent editing (default 0.025).
#' @param theta_zero_p probability that a cell type does not express a site's host gene (default 0.25).
#' @param coverage `"realistic"` (default) or `"poisson"`.
#' @param mean_coverage median (realistic) or mean (poisson) coverage (default 40).
#' @param seed random seed (default 0).
#' @return a list with `bulk_editing`, `coverage`, `bulk_expression` (sites x donors), `genotypes` (variants x
#'   donors), `proportions` (donors x cell types), `theta` (sites x cell types), `pairs` (data.frame site_id,
#'   variant_id) and `truth`: `site_type` (named by site), `mu`, `sd`, `beta` (sites x cell types; beta is the effect of
#'   the causal variant), `latent` (named list per cell type of sites x donors true editing, NA where the cell type
#'   does not express the host gene), `causal_variant` (named by site) and `variants` (data.frame per variant: site_id,
#'   variant_id, role, r2 with the causal variant).
#' @examples
#' co <- simulate_edqtl_cohort(n_donors = 300, n_sites_per_type = 2, n_variants = 3)
#' dim(co$bulk_editing); dim(co$genotypes); head(co$pairs)
#' @export
simulate_edqtl_cohort <- function(n_donors = 1000, celltypes = c("Neutrophils", "Monocytes", "CD4", "NK"),
                                  composition = c(0.45, 0.23, 0.18, 0.14), concentration = 20, n_sites_per_type = 10,
                                  n_variants = 5, maf = 0.25, beta = 0.05, baseline = c(0.05, 0.30), sd = 0.025,
                                  theta_zero_p = 0.25, coverage = c("realistic", "poisson"), mean_coverage = 40,
                                  seed = 0) {
  coverage <- match.arg(coverage)
  K <- length(celltypes)
  if (length(composition) != K) stop("composition must have one value per cell type", call. = FALSE)
  if (n_variants < 1) stop("n_variants must be >= 1", call. = FALSE)
  set.seed(seed)
  donors <- sprintf("donor%04d", seq_len(n_donors))
  P <- matrix(stats::rgamma(n_donors * K, shape = rep(composition * concentration, n_donors)), n_donors, K, byrow = TRUE)
  P <- P / rowSums(P); dimnames(P) <- list(donors, celltypes)
  types <- rep(c(paste0("only_", celltypes), "shared", "null"), each = n_sites_per_type)
  S <- length(types)
  # unique, realistic-looking site ids (chrom:pos:strand)
  chrom <- sample(as.character(1:22), S, replace = TRUE); pos <- sample.int(2.4e8, S); strand <- sample(c("+", "-"), S, replace = TRUE)
  sites <- format_site_id(chrom, pos, strand)
  depth <- if (coverage == "realistic") exp(stats::rnorm(n_donors, 0, 0.6)) else NULL
  var_ids <- as.vector(t(outer(sites, seq_len(n_variants), function(s, j) paste0(s, "_v", j))))
  G <- matrix(NA_real_, length(var_ids), n_donors, dimnames = list(var_ids, donors))
  Y <- CV <- matrix(NA_real_, S, n_donors, dimnames = list(sites, donors))
  TH <- MU <- SD <- BE <- matrix(0, S, K, dimnames = list(sites, celltypes))
  LAT <- stats::setNames(lapply(celltypes, function(h) matrix(NA_real_, S, n_donors, dimnames = list(sites, donors))), celltypes)
  vrole <- character(length(var_ids)); vr2 <- numeric(length(var_ids))
  for (s in seq_len(S)) {
    ty <- types[s]
    aff <- switch(ty, shared = celltypes, null = character(0), sub("^only_", "", ty))
    th <- exp(stats::rnorm(K, 0, 0.5))
    if (theta_zero_p > 0) {
      off <- stats::runif(K) < theta_zero_p; off[celltypes %in% aff & ty != "shared"] <- FALSE
      while (sum(!off) < 2) off[sample(which(off), 1)] <- FALSE
      th[off] <- 0
    }
    b <- ifelse(celltypes %in% aff & th > 0, beta, 0)
    mu <- stats::runif(1, baseline[1], baseline[2]) * stats::runif(K, 0.6, 1.4)
    sds <- sd * stats::runif(K, 0.5, 1.5)
    # variants: causal from two haplotypes, two in LD, the rest independent
    h1 <- stats::rbinom(n_donors, 1, maf); h2 <- stats::rbinom(n_donors, 1, maf)
    gs <- matrix(0, n_variants, n_donors); role <- character(n_variants)
    for (j in seq_len(n_variants)) {
      if (j == 1) { gs[j, ] <- h1 + h2; role[j] <- "causal" }
      else if (j <= 3) {
        cp <- c(0.9, 0.6)[j - 1]
        a1 <- ifelse(stats::runif(n_donors) < cp, h1, stats::rbinom(n_donors, 1, maf))
        a2 <- ifelse(stats::runif(n_donors) < cp, h2, stats::rbinom(n_donors, 1, maf))
        gs[j, ] <- a1 + a2; role[j] <- paste0("ld", j - 1)
      } else { gs[j, ] <- stats::rbinom(n_donors, 2, stats::runif(1, 0.05, 0.5)); role[j] <- "independent" }
    }
    g <- gs[1, ]
    Z <- pmin(pmax(vapply(seq_len(K), function(h) mu[h] + b[h] * g + stats::rnorm(n_donors, 0, sds[h]), numeric(n_donors)), 0), 1)
    W <- sweep(P, 2, th, `*`); phi <- W / rowSums(W)
    cv <- if (coverage == "realistic") {
      sdp <- min(max(exp(stats::rnorm(1, log(mean_coverage), 1)), 5), 1000)
      pmax(stats::rpois(n_donors, sdp * depth * rowSums(W) / mean(rowSums(W))), 2)
    } else stats::rpois(n_donors, mean_coverage) + 10
    Y[s, ] <- stats::rbinom(n_donors, cv, pmin(pmax(rowSums(phi * Z), 0), 1)) / cv; CV[s, ] <- cv
    vs <- (s - 1) * n_variants + seq_len(n_variants); G[vs, ] <- gs; vrole[vs] <- role
    vr2[vs] <- apply(gs, 1, function(x) if (stats::var(x) > 0) stats::cor(x, g)^2 else NA_real_)
    TH[s, ] <- th; MU[s, ] <- mu; SD[s, ] <- sds; BE[s, ] <- b
    for (h in seq_len(K)) if (th[h] > 0) LAT[[h]][s, ] <- Z[, h]
  }
  list(bulk_editing = Y, coverage = CV, bulk_expression = simulate_bulk_expression(P, TH, seed = seed + 1),
       genotypes = G, proportions = P, theta = TH,
       pairs = data.frame(site_id = rep(sites, each = n_variants), variant_id = var_ids, stringsAsFactors = FALSE),
       truth = list(site_type = stats::setNames(types, sites), mu = MU, sd = SD, beta = BE, latent = LAT,
                    causal_variant = stats::setNames(paste0(sites, "_v1"), sites),
                    variants = data.frame(site_id = rep(sites, each = n_variants), variant_id = var_ids, role = vrole,
                                          r2_with_causal = vr2, stringsAsFactors = FALSE)))
}
