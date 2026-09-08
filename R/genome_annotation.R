#' Format a standardized caEditR site id
#'
#' Every function in this package that takes or returns RNA-editing site
#' ids (`bulk_editing`/`coverage` row names, `map_sites_to_genes()`,
#' `build_coverage_from_expression()`, `simulate_reference_and_cohort()`,
#' etc.) uses exactly ONE standardized id format: `"chrom:pos"` or
#' `"chrom:pos:strand"` (e.g. `"10:100232436:-"`) -- chromosome name
#' WITHOUT a leading "chr" (bare Ensembl/GENCODE-style, matching this
#' project's own real editing-call output), a 1-based genomic position,
#' and an optional strand (`"+"`/`"-"`; omit or use `"*"`/`NA` if unknown
#' or not strand-specific). This is the canonical way to construct one --
#' use it instead of hand-building the string yourself, so the format
#' always matches exactly what `map_sites_to_genes()`/
#' `build_coverage_from_expression()` (and everything else in this
#' package) expect. (`map_sites_to_genes()` also ACCEPTS "chr"-prefixed
#' chromosome names on input, auto-detected, for convenience -- but this
#' function always OUTPUTS the bare, no-"chr" canonical form.)
#'
#' @param chrom chromosome name(s), with or without a leading "chr" (e.g. `"10"` or `"chr10"`).
#' @param pos integer genomic position(s), 1-based.
#' @param strand optional strand(s), `"+"`/`"-"`/`"*"`/`NA` (default `"*"` = unknown/unspecified).
#' @return character vector of standardized site ids, recycled to the
#'   longest of `chrom`/`pos`/`strand`.
#' @examples
#' format_site_id("10", 100232436, "-")
#' format_site_id("chr10", c(100232436, 100232925), c("-", "-"))
#' format_site_id("1", 12831014)  # strand unspecified -> "1:12831014", no trailing ":*"
#' @export
format_site_id <- function(chrom, pos, strand = "*") {
  chrom <- sub("^chr", "", chrom)
  n <- max(length(chrom), length(pos), length(strand))
  chrom <- rep_len(chrom, n)
  pos <- rep_len(pos, n)
  strand <- rep_len(strand, n)
  if (!all(strand %in% c("+", "-", "*") | is.na(strand))) {
    stop('strand must be "+", "-", "*", or NA', call. = FALSE)
  }
  ifelse(is.na(strand) | strand == "*",
         sprintf("%s:%d", chrom, as.integer(pos)),
         sprintf("%s:%d:%s", chrom, as.integer(pos), strand))
}

#' Parse this package's own "chrom:pos" / "chrom:pos:strand" site ids
#' (see `format_site_id()` for the canonical way to construct one).
#' @keywords internal
.parse_site_ids <- function(site_ids) {
  parts <- strsplit(site_ids, ":", fixed = TRUE)
  if (any(lengths(parts) < 2)) {
    bad <- site_ids[lengths(parts) < 2][1]
    stop("Every site id must be at least 'chrom:pos' (optionally 'chrom:pos:strand'); ",
         "found a malformed id: '", bad, "'", call. = FALSE)
  }
  chrom <- vapply(parts, `[`, character(1), 1)
  pos <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 2)))
  if (anyNA(pos)) stop("Every site id's position must be an integer (chrom:pos[:strand]).", call. = FALSE)
  strand <- vapply(parts, function(p) if (length(p) >= 3 && p[3] %in% c("+", "-")) p[3] else "*", character(1))
  data.frame(site_id = site_ids, chrom = chrom, pos = pos, strand = strand, stringsAsFactors = FALSE)
}

#' Which UCSC `TxDb` annotation package covers a given genome build.
#' @keywords internal
.txdb_pkg_for_genome <- function(genome) {
  genome <- match.arg(genome, c("hg19", "hg38"))
  switch(genome,
    hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene"
  )
}

#' Ensure the genome-annotation packages needed by `map_sites_to_genes()`
#' are installed for a given build (`hg19` or `hg38`).
#'
#' Routes through the one shared `.ensure_bioc_installed()` helper
#' (`ensure_bioc.R`) -- same single-central-installer design as
#' `.ensure_music_stack_installed()`, not a separate reimplementation.
#' @return the `TxDb` package name for `genome` (invisibly usable to `get()` the TxDb object).
#' @keywords internal
.ensure_genome_annotation_installed <- function(genome) {
  txdb_pkg <- .txdb_pkg_for_genome(genome)
  .ensure_bioc_installed(c("GenomicFeatures", "GenomicRanges", "IRanges", "S4Vectors",
                            "AnnotationDbi", "org.Hs.eg.db", txdb_pkg))
  txdb_pkg
}

#' Map RNA-editing site ids to their host gene (Ensembl id), on hg19 or hg38
#'
#' Uses the real Bioconductor UCSC `TxDb` gene-model annotation for the
#' requested genome build (`TxDb.Hsapiens.UCSC.hg19.knownGene` or
#' `..hg38..`), plus `org.Hs.eg.db` to translate its Entrez gene ids to
#' Ensembl gene ids -- not a new/custom gene-annotation source. Both
#' annotation packages (plus their `GenomicFeatures`/`GenomicRanges`
#' dependencies) are installed automatically on first use if missing (see
#' `.ensure_genome_annotation_installed()`), same auto-install design as
#' `TCA_Like()`/`estimate_proportions_music()` elsewhere in this package.
#'
#' Each site is matched to any gene whose body overlaps its genomic
#' position. If a site's own strand is known (`"chrom:pos:strand"`, this
#' package's own site-id convention) and more than one gene overlaps at
#' that position, a gene on the SAME strand is preferred (common for
#' ADAR-edited sites in overlapping sense/antisense transcripts); ties are
#' broken by taking the first match, deterministically.
#'
#' @param site_ids character vector of standardized caEditR site ids --
#'   `"chrom:pos"` or `"chrom:pos:strand"` (e.g. `"10:100232436:-"`; see
#'   `format_site_id()` for the canonical way to construct these).
#'   Chromosome names with or without a leading "chr" are both accepted.
#' @param genome `"hg19"` or `"hg38"` -- which genome build `site_ids`'
#'   coordinates are on (REQUIRED -- the same numeric position means a
#'   different gene, or no gene at all, on each build; see the
#'   `map_sites_to_genes()` unit tests for a direct example of this).
#' @return a data.frame, one row per site id, columns `site_id`, `chrom`,
#'   `pos`, `strand`, `entrez_gene_id`, `ensembl_gene_id` (the last two are
#'   `NA` for a site with no overlapping gene).
#' @examples
#' # Real site ids bundled with this package (GRCh38/hg38 coordinates --
#' # see this project's own README for the alignment pipeline's genome build).
#' extdata <- system.file("extdata", package = "caEditR")
#' bulk <- read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1)
#' map_sites_to_genes(rownames(bulk)[1:3], genome = "hg38")
#' @export
map_sites_to_genes <- function(site_ids, genome = c("hg19", "hg38")) {
  genome <- match.arg(genome)
  txdb_pkg <- .ensure_genome_annotation_installed(genome)
  txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))

  parsed <- .parse_site_ids(site_ids)
  ucsc_chrom <- ifelse(startsWith(parsed$chrom, "chr"), parsed$chrom, paste0("chr", parsed$chrom))

  site_gr <- GenomicRanges::GRanges(
    seqnames = ucsc_chrom,
    ranges = IRanges::IRanges(start = parsed$pos, width = 1),
    strand = parsed$strand
  )
  # suppressMessages(): genes() reports how many genes it drops for having
  # exons on multiple strands/seqnames (can't be a single GRange) -- a
  # benign, expected note about GenomicFeatures' own representation limits,
  # not something caEditR needs to surface every call.
  gene_gr <- suppressMessages(GenomicFeatures::genes(txdb))  # names(gene_gr) == Entrez gene ids

  hits <- GenomicRanges::findOverlaps(site_gr, gene_gr, ignore.strand = TRUE)
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)

  entrez_ids <- rep(NA_character_, length(site_ids))
  if (length(q) > 0) {
    site_strand <- as.character(GenomicRanges::strand(site_gr))[q]
    gene_strand <- as.character(GenomicRanges::strand(gene_gr))[s]
    strand_match <- site_strand == "*" | site_strand == gene_strand

    for (idx in split(seq_along(q), q)) {
      keep <- if (any(strand_match[idx])) idx[strand_match[idx]][1] else idx[1]
      entrez_ids[q[keep]] <- names(gene_gr)[s[keep]]
    }
  }

  unique_entrez <- unique(stats::na.omit(entrez_ids))
  ensembl_lookup <- stats::setNames(character(0), character(0))
  if (length(unique_entrez) > 0) {
    org_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
    mapping <- suppressMessages(AnnotationDbi::select(
      org_db, keys = unique_entrez, keytype = "ENTREZID", columns = "ENSEMBL"
    ))
    mapping <- mapping[!is.na(mapping$ENSEMBL) & !duplicated(mapping$ENTREZID), ]
    ensembl_lookup <- stats::setNames(mapping$ENSEMBL, mapping$ENTREZID)
  }

  data.frame(
    site_id = site_ids, chrom = parsed$chrom, pos = parsed$pos, strand = parsed$strand,
    entrez_gene_id = entrez_ids,
    ensembl_gene_id = unname(ensembl_lookup[entrez_ids]),
    stringsAsFactors = FALSE
  )
}

#' Resolve the `coverage` argument shared by `caRD_edit()`/`caNRD_edit()`:
#' either the caller supplied a real coverage matrix directly (preferred,
#' returned as-is), or supplied `expression` instead, in which case
#' `build_coverage_from_expression()` derives one -- the ONE place this
#' fallback logic lives, used by both functions, not duplicated per caller.
#' @keywords internal
.resolve_coverage <- function(bulk_editing, coverage, expression, genome, coverage_scale, unmapped_floor) {
  if (is.null(coverage) && is.null(expression)) {
    stop("Must supply either `coverage` (a real sites x samples read-depth ",
         "matrix) or `expression` (a genes x samples matrix, to derive ",
         "coverage from automatically via build_coverage_from_expression()).",
         call. = FALSE)
  }
  if (!is.null(coverage)) {
    if (!is.null(expression)) {
      message("Both `coverage` and `expression` were supplied -- using `coverage` ",
              "directly (real coverage always takes precedence); `expression` is ignored.")
    }
    return(as.matrix(coverage))
  }
  genome <- match.arg(genome, c("hg19", "hg38"))
  build_coverage_from_expression(rownames(bulk_editing), expression, genome = genome,
                                  coverage_scale = coverage_scale, unmapped_floor = unmapped_floor)
}

#' Build a per-site coverage matrix from per-sample gene expression
#'
#' Real RNA-editing read coverage at a site scales with how highly its
#' host gene is expressed in that sample -- this approximates a per-site,
#' per-sample coverage matrix directly from an existing per-sample gene
#' expression matrix (Ensembl gene ids), via `map_sites_to_genes()`, for
#' use as the `coverage` argument to `caRD_edit()`/`caNRD_edit()` when a
#' real per-site read-coverage matrix (e.g. from an editing caller) isn't
#' available. This is an approximation, not real read coverage -- prefer
#' real per-site coverage from your editing calls whenever you have it.
#'
#' @param site_ids character vector of standardized caEditR site ids --
#'   `"chrom:pos"` or `"chrom:pos:strand"` (see `format_site_id()`) -- the
#'   rows of the coverage matrix to build.
#' @param expression numeric matrix, genes (rows, Ensembl gene ids --
#'   version suffixes like `"ENSG00000141510.16"` are stripped
#'   automatically before matching) x samples (columns).
#' @param genome `"hg19"` or `"hg38"` -- which genome build `site_ids`'
#'   coordinates are on (REQUIRED -- see `map_sites_to_genes()`).
#' @param coverage_scale multiply each site's host-gene expression by this
#'   factor before using it as coverage (default 1 -- set based on your own
#'   expression units; e.g. if `expression` is in TPM and you want
#'   read-count-like magnitudes, scale up accordingly).
#' @param unmapped_floor coverage assigned to sites whose gene couldn't be
#'   determined, or whose gene isn't present in `expression` (default 1 --
#'   deliberately low so `caRD_edit()`/`caNRD_edit()`'s own `min_coverage`
#'   threshold correctly flags these as low-confidence rather than being
#'   silently trusted).
#' @return numeric matrix, sites (rows, `site_ids`) x samples (columns,
#'   `colnames(expression)`), with a `"site_to_gene"` attribute (the full
#'   `map_sites_to_genes()` result, for inspection/diagnostics).
#' @examples
#' extdata <- system.file("extdata", package = "caEditR")
#' bulk <- read.csv(file.path(extdata, "example_bulk_editing_ratios.csv"), row.names = 1)
#' site_ids <- rownames(bulk)[1:3]
#' # A toy expression matrix keyed by the actual genes those 3 sites fall
#' # in (see map_sites_to_genes(site_ids, "hg38")$ensembl_gene_id) --
#' # substitute your own real per-sample Ensembl-indexed expression matrix.
#' genes <- map_sites_to_genes(site_ids, genome = "hg38")$ensembl_gene_id
#' expression <- matrix(c(50, 200, 30, 400, 20, 5), nrow = length(genes),
#'                       dimnames = list(genes, c("sampleA", "sampleB")))
#' build_coverage_from_expression(site_ids, expression, genome = "hg38")
#' @export
build_coverage_from_expression <- function(site_ids, expression, genome = c("hg19", "hg38"),
                                            coverage_scale = 1, unmapped_floor = 1) {
  genome <- match.arg(genome)
  expression <- as.matrix(expression)
  sample_ids <- colnames(expression)
  if (is.null(sample_ids)) stop("expression must have column names (sample ids)", call. = FALSE)

  site_to_gene <- map_sites_to_genes(site_ids, genome = genome)
  expr_gene_ids <- sub("\\..*$", "", rownames(expression))  # strip Ensembl version suffixes
  match_row <- match(site_to_gene$ensembl_gene_id, expr_gene_ids)

  coverage <- matrix(unmapped_floor, nrow = length(site_ids), ncol = length(sample_ids),
                      dimnames = list(site_ids, sample_ids))
  mapped <- !is.na(match_row)
  if (any(mapped)) {
    coverage[mapped, ] <- pmax(round(expression[match_row[mapped], , drop = FALSE] * coverage_scale), unmapped_floor)
  }

  message(sprintf(
    "build_coverage_from_expression(): %d/%d sites mapped to a gene present in `expression` (genome=%s); %d unmapped site(s) got the floor coverage (%g).",
    sum(mapped), length(site_ids), genome, length(site_ids) - sum(mapped), unmapped_floor
  ))

  attr(coverage, "site_to_gene") <- site_to_gene
  coverage
}
