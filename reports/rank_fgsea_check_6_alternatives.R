# Notebook 05, Exercise 1, Step 5 - cutoff-based test and Reactome, side by side

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fgsea)
  library(msigdbr)
})

genes_by_outcome <- read_csv(here("data", "genes_by_outcome.csv"), show_col_types = FALSE)

stat_vec <- genes_by_outcome$stat
names(stat_vec) <- genes_by_outcome$feature
stat_vec <- sort(stat_vec[!is.na(stat_vec)], decreasing = TRUE)
universe <- names(stat_vec)
N <- length(universe)

# msigdbr, built the safe way: distinct gene-set pairs before splitting
h_all <- if (packageVersion("msigdbr") >= "10.0.0") {
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
} else {
  msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
}
reactome <- h_all |>
  distinct(gs_name, ensembl_gene) |>
  with(split(ensembl_gene, gs_name))

h_hm <- if (packageVersion("msigdbr") >= "10.0.0") {
  msigdbr(species = "Homo sapiens", collection = "H")
} else {
  msigdbr(species = "Homo sapiens", category = "H")
}
hallmark <- h_hm |>
  distinct(gs_name, ensembl_gene) |>
  with(split(ensembl_gene, gs_name))

# keep only genes we actually tested, apply the same 15-500 filter everywhere
in_universe <- function(sets) {
  sets |>
    map(~ intersect(.x, universe)) |>
    keep(~ length(.x) >= 15 & length(.x) <= 500)
}
hallmark_u  <- in_universe(hallmark)
reactome_u  <- in_universe(reactome)

# ---- method A: hypergeometric on significant genes, up and down separately ----
# "Of the significant genes, is this set over-represented among the up ones?"
hyper <- function(sets, sig_ids) {
  n <- length(sig_ids)
  map_dfr(names(sets), function(nm) {
    m  <- length(sets[[nm]])            # set genes in the universe
    k  <- sum(sets[[nm]] %in% sig_ids)  # overlap with the significant list
    p  <- phyper(k - 1, m, N - m, n, lower.tail = FALSE)
    tibble(pathway = nm, overlap = k, set_size = m, n_sig = n, p = p)
  }) |>
    mutate(padj = p.adjust(p, method = "BH"))
}

sig_up   <- genes_by_outcome$feature[genes_by_outcome$padj < 0.05 & genes_by_outcome$stat > 0]
sig_down <- genes_by_outcome$feature[genes_by_outcome$padj < 0.05 & genes_by_outcome$stat < 0]

hyp_up   <- hyper(hallmark_u, sig_up)
hyp_down <- hyper(hallmark_u, sig_down)

cat("== Hypergeometric (cutoff-based), Hallmark ==\n")
cat("significant in UP genes (padj<0.05):", sum(hyp_up$padj < 0.05), "\n")
cat("significant in DOWN genes (padj<0.05):", sum(hyp_down$padj < 0.05), "\n\n")
cat("Top 5 up:\n");   print(hyp_up   |> arrange(padj) |> select(pathway, overlap, set_size, padj) |> head(5))
cat("\nTop 5 down:\n"); print(hyp_down |> arrange(padj) |> select(pathway, overlap, set_size, padj) |> head(5))

# ---- method B: fgsea on Reactome ----------------------------------------------
set.seed(42)
fg_react <- fgsea(stats = stat_vec, pathways = reactome_u, minSize = 15, maxSize = 500)
fg_react <- as.data.frame(fg_react) |> arrange(padj) |> as_tibble()

cat("\n== fgsea on Reactome ==\n")
cat("sets tested:", nrow(fg_react),
    "| significant:", sum(fg_react$padj < 0.05),
    "| higher in preterm:", sum(fg_react$padj < 0.05 & fg_react$NES > 0), "\n")
print(fg_react |> select(pathway, size, NES, padj) |> head(10) |> mutate(across(where(is.numeric), \(x) signif(x, 3))))

# ---- side by side with the Hallmark fgsea result ------------------------------
fg_hall <- read_csv(here("data", "fgsea_hallmark_outcome.csv"), show_col_types = FALSE)

compare <- fg_hall |>
  select(pathway, fgsea_NES = NES, fgsea_padj = padj) |>
  full_join(hyp_up   |> transmute(pathway, hyper_up_padj = padj, hyper_up_overlap = overlap), by = "pathway") |>
  full_join(hyp_down |> transmute(pathway, hyper_down_padj = padj, hyper_down_overlap = overlap), by = "pathway") |>
  mutate(
    fgsea_sig        = fgsea_padj < 0.05,
    hyper_up_sig     = hyper_up_padj < 0.05,
    hyper_down_sig   = hyper_down_padj < 0.05,
    direction_agrees = case_when(
      fgsea_sig & hyper_up_sig & fgsea_NES > 0      ~ "agree: preterm",
      fgsea_sig & hyper_down_sig & fgsea_NES < 0    ~ "agree: term",
      fgsea_sig & !hyper_up_sig & !hyper_down_sig    ~ "fgsea only",
      !fgsea_sig & (hyper_up_sig | hyper_down_sig)   ~ "hypergeom only",
      TRUE                                           ~ "neither"
    )
  ) |>
  arrange(desc(fgsea_NES))

cat("\n== Hallmark: fgsea vs hypergeometric, all 50 sets ==\n")
print(compare |> count(direction_agrees) |> arrange(n))

cat("\nSets significant by BOTH methods:\n")
print(compare |> filter(fgsea_sig, hyper_up_sig | hyper_down_sig) |>
      select(pathway, fgsea_NES, fgsea_padj, hyper_up_padj, hyper_down_padj, direction_agrees))

write_csv(compare, here("data", "hallmark_methods_side_by_side.csv"))
write_csv(fg_react |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")),
          here("data", "fgsea_reactome_outcome.csv"))
cat("\nSaved: data/hallmark_methods_side_by_side.csv, data/fgsea_reactome_outcome.csv\n")
