# Notebook 04, Exercise 2, Step 5 - the ~ outcome test on genes
# Produces data/genes_by_outcome.csv, which notebook 05's setup chunk reads.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
})

source(here("R/count_tables.R"))
d <- read_counts("instructional")

# --- week 18 samples only -----------------------------------------------------
meta18 <- d$metadata |> filter(time_point == "week_18")
counts18 <- d$genes[, meta18$sample_id, drop = FALSE]

# --- the order check: count columns must match metadata rows exactly ----------
stopifnot(identical(colnames(counts18), meta18$sample_id))
cat("column order check passed:", ncol(counts18), "samples\n")

# --- gene filter: at least 10 reads in at least a quarter of the samples ------
keep <- rowSums(counts18 >= 10) >= ceiling(nrow(meta18) / 4)
genes_keep <- counts18[keep, , drop = FALSE]
cat("genes after filter:", nrow(genes_keep), "\n")

# --- DESeq2 on RAW counts -----------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(genes_keep),
  colData   = as.data.frame(meta18),
  design    = ~ outcome
)
dds$outcome <- relevel(factor(dds$outcome), ref = "term")  # contrast: preterm vs term
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("outcome", "preterm", "term")))
res$feature <- rownames(res)

# --- join symbols and gene names ---------------------------------------------
out <- res |>
  as_tibble() |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, symbol, gene_name) |>
  arrange(padj)

dir.create(here("data"), showWarnings = FALSE)
write_csv(out, here("data", "genes_by_outcome.csv"))

cat("genes tested:", sum(!is.na(out$pvalue)), "\n")
cat("padj < 0.05 :", sum(out$padj < 0.05, na.rm = TRUE), "\n")
cat("positive stat (higher in preterm):", sum(out$stat > 0, na.rm = TRUE), "\n")
