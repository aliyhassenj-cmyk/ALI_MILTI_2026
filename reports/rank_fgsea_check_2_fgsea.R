# Notebook 05, Exercise 1, Step 1 - rank genes by stat, fgsea vs Hallmark
# Plus the identifier check the user asked for.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fgsea)
  library(msigdbr)
})

genes_by_outcome <- read_csv(here("data", "genes_by_outcome.csv"), show_col_types = FALSE)

# ---- the ranking vector: stat, named by Ensembl ID, NAs dropped --------------
stat_vec <- genes_by_outcome$stat
names(stat_vec) <- genes_by_outcome$feature
stat_vec <- stat_vec[!is.na(stat_vec)]
stat_vec <- sort(stat_vec, decreasing = TRUE)
cat("genes in ranking:", length(stat_vec), "\n")

# ---- Hallmark sets from msigdbr, Ensembl IDs --------------------------------
if (packageVersion("msigdbr") >= "10.0.0") {
  h <- msigdbr(species = "Homo sapiens", collection = "H")
} else {
  h <- msigdbr(species = "Homo sapiens", category = "H")
}
cat("msigdbr columns:", paste(names(h), collapse = ", "), "\n")
ens_col <- intersect(c("gs_ensembl_gene", "ensembl_gene"), names(h))[1]
cat("using Ensembl column:", ens_col, "\n")

hallmark <- split(h[[ens_col]], h$gs_name)
sets <- lapply(hallmark, unique)

# ---- fgsea -------------------------------------------------------------------
set.seed(42)
fg <- fgsea(stats = stat_vec, pathways = sets, minSize = 15, maxSize = 500)
fg <- as.data.frame(fg) |> arrange(padj) |> as_tibble()

cat("\n== fgsea results ==\n")
cat("sets tested:", nrow(fg), "\n")
cat("significant at padj < 0.05:", sum(fg$padj < 0.05), "\n")
cat("  higher in PRETERM (NES > 0):", sum(fg$padj < 0.05 & fg$NES > 0), "\n")
cat("  higher in TERM    (NES < 0):", sum(fg$padj < 0.05 & fg$NES < 0), "\n")

cat("\nTop 10 by padj:\n")
print(fg |> select(pathway, size, NES, padj, leadingEdge) |>
  mutate(leadingEdge = map_int(leadingEdge, length)) |>
  head(10))

write_csv(fg |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")),
          here("data", "fgsea_hallmark_outcome.csv"))
