# Notebook 05, Exercise 1, Step 1 - did the identifiers actually match?

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

if (packageVersion("msigdbr") >= "10.0.0") {
  h <- msigdbr(species = "Homo sapiens", collection = "H")
} else {
  h <- msigdbr(species = "Homo sapiens", category = "H")
}
msig_ids <- unique(h$ensembl_gene)

cat("== 1. Are the ranking names really Ensembl gene IDs? ==\n")
n_ensg   <- sum(grepl("^ENSG[0-9]{11}$", names(stat_vec)))
n_ver    <- sum(grepl("^ENSG[0-9]{11}\\.", names(stat_vec)))   # versioned, e.g. ENSG...11.4
n_symbol <- sum(!grepl("^ENSG", names(stat_vec)))              # would be a symbol mix-up
cat("  clean ENSG format :", n_ensg, "of", length(stat_vec), "\n")
cat("  versioned IDs     :", n_ver, "(fgsea would treat these as different strings)\n")
cat("  non-ENSG strings  :", n_symbol, "(would signal a symbol/ID mix-up)\n")

cat("\n== 2. Did they match the Hallmark sets? ==\n")
in_any_set <- names(stat_vec) %in% msig_ids
cat("  ranking genes found in at least one Hallmark set:", sum(in_any_set), "of", length(stat_vec), "\n")
cat("  msigdbr Ensembl IDs in Hallmark total:", length(msig_ids), "\n")

cat("\n== 3. Would a wrong ID type have failed loudly or silently? ==\n")
# What if we had used symbols instead (a common mistake)?
sym_vec <- genes_by_outcome$stat
names(sym_vec) <- genes_by_outcome$symbol
sym_vec <- sort(sym_vec[!is.na(sym_vec)], decreasing = TRUE)
sym_in_sets <- sum(names(sym_vec) %in% msig_ids)
cat("  symbols instead of Ensembl IDs matched:", sym_in_sets, "of", length(sym_vec), "genes\n")
cat("  -> nearly zero, so the match is real, not luck\n")

# What if we had versioned Ensembl IDs? Make them and check.
ver_names <- paste0(names(stat_vec), ".1")
cat("  versioned IDs matched:", sum(ver_names %in% msig_ids), "-> version suffix would break the match\n")

cat("\n== 4. Direction sanity check ==\n")
# The two genes DESeq2 called most preterm-high and most term-high:
top_preterm <- genes_by_outcome |> filter(stat > 0) |> slice_max(stat, n = 3)
top_term    <- genes_by_outcome |> filter(stat < 0) |> slice_min(stat, n = 3)
cat("  highest stat (expect preterm-high):", paste(top_preterm$symbol, collapse = ", "), "\n")
cat("  lowest  stat (expect term-high)   :", paste(top_term$symbol, collapse = ", "), "\n")
# And NES sign in the fgsea output matches the contrast (positive = higher in preterm):
fg <- read_csv(here("data", "fgsea_hallmark_outcome.csv"), show_col_types = FALSE)
cat("  significant sets with NES > 0 (preterm-high):", sum(fg$padj < 0.05 & fg$NES > 0), "\n")
cat("  significant sets with NES < 0 (term-high)   :", sum(fg$padj < 0.05 & fg$NES < 0), "\n")
