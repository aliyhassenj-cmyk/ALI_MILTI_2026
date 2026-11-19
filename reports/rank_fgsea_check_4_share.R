# Notebook 05, Exercise 1, Step 3 - interpret the top sets, then shared genes

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(msigdbr)
})

fg <- read_csv(here("data", "fgsea_hallmark_outcome.csv"), show_col_types = FALSE) |>
  mutate(leadingEdge = str_split(leadingEdge, ";"))

top5 <- fg |> arrange(padj) |> slice_head(n = 5)

if (packageVersion("msigdbr") >= "10.0.0") {
  h <- msigdbr(species = "Homo sapiens", collection = "H")
} else {
  h <- msigdbr(species = "Homo sapiens", category = "H")
}
# msigdbr has several rows per gene-set pair (different database mappings),
# so take distinct pairs BEFORE splitting, or sizes come out wrong.
sets <- h |> distinct(gs_name, ensembl_gene) |> with(split(ensembl_gene, gs_name))

cat("== Top five sets (by padj) ==\n")
print(top5 |> select(pathway, size, NES, padj) |> mutate(across(where(is.numeric), signif, 3)))

cat("\n== Pairwise shared genes between the top five sets ==\n")
pairs <- utils::combn(top5$pathway, 2)
pair_tab <- map_dfr(seq_len(ncol(pairs)), function(i) {
  a <- pairs[1, i]; b <- pairs[2, i]
  sa <- sets[[a]]; sb <- sets[[b]]
  shared <- intersect(sa, sb)
  tibble(
    set_a    = str_remove(a, "HALLMARK_"),
    set_b    = str_remove(b, "HALLMARK_"),
    size_a   = length(sa),
    size_b   = length(sb),
    shared   = length(shared),
    pct_of_a = round(100 * shared / length(sa)),
    pct_of_b = round(100 * shared / length(sb)),
    jaccard  = round(shared / length(union(sa, sb)), 2)
  )
}) |> arrange(desc(jaccard))

print(pair_tab, n = 10, width = 120)

cat("\n== Shared genes between LEADING EDGES (the genes doing the work) ==\n")
le <- set_names(top5$leadingEdge, str_remove(top5$pathway, "HALLMARK_"))
le_pairs <- utils::combn(names(le), 2)
for (i in seq_len(ncol(le_pairs))) {
  a <- le_pairs[1, i]; b <- le_pairs[2, i]
  sh <- intersect(le[[a]], le[[b]])
  cat(sprintf("%-22s vs %-22s: %3d of %3d/%3d leading-edge genes shared\n",
              a, b, length(sh), length(le[[a]]), length(le[[b]])))
}

all_le <- unlist(le)
shared_count <- table(all_le)
cat("\nGenes in the most leading edges at once:\n")
print(sort(shared_count[shared_count >= 3], decreasing = TRUE))
