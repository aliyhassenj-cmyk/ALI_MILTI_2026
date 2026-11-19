# Notebook 04, Exercise 2, Steps 1+3 - species ~ outcome, poscounts version
# Produces data/species_by_outcome.csv, which notebook 05 Exercise 2 reads.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
})

source(here("R/count_tables.R"))
d <- read_counts("instructional")

meta18 <- d$metadata |> filter(time_point == "week_18")
sp18 <- d$species[, meta18$sample_id, drop = FALSE]
stopifnot(identical(colnames(sp18), meta18$sample_id))   # order check

# species present in at least 10% of samples
keep <- rowSums(sp18 > 0) >= 0.10 * nrow(meta18)
sp_keep <- sp18[keep, , drop = FALSE]
cat("species after 10% presence filter:", nrow(sp_keep), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(sp_keep),
  colData   = as.data.frame(meta18),
  design    = ~ outcome
)
dds$outcome <- relevel(factor(dds$outcome), ref = "term")
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

res <- as.data.frame(results(dds, contrast = c("outcome", "preterm", "term")))
res$feature <- rownames(res)

out <- res |>
  as_tibble() |>
  left_join(select(d$taxonomy, taxon, Genus, Family),
            by = c("feature" = "taxon")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, Genus, Family) |>
  arrange(padj)

write_csv(out, here("data", "species_by_outcome.csv"))
cat("species tested:", sum(!is.na(out$pvalue)), "\n")
cat("padj < 0.05   :", sum(out$padj < 0.05, na.rm = TRUE), "\n")
print(head(out, 5))
