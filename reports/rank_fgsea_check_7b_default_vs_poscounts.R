# Compare default vs poscounts size factors for the species ~ outcome test
# (notebook 04, Exercise 2, Steps 1 vs 3, side by side)

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(DESeq2)
})

source(here("R/count_tables.R"))
d <- read_counts("instructional")

meta18 <- d$metadata |> filter(time_point == "week_18")
sp18   <- d$species[, meta18$sample_id, drop = FALSE]
stopifnot(identical(colnames(sp18), meta18$sample_id))

keep <- rowSums(sp18 > 0) >= 0.10 * nrow(meta18)
sp_keep <- sp18[keep, , drop = FALSE]
cat("species after filter:", nrow(sp_keep), "\n")

run <- function(type = c("default", "poscounts")) {
  type <- match.arg(type)
  dds <- DESeqDataSetFromMatrix(round(sp_keep), as.data.frame(meta18), ~ outcome)
  dds$outcome <- relevel(factor(dds$outcome), ref = "term")
  if (type == "poscounts") {
    dds <- estimateSizeFactors(dds, type = "poscounts")
  }
  dds <- DESeq(dds)
  r <- as.data.frame(results(dds, contrast = c("outcome", "preterm", "term")))
  r$feature <- rownames(r)
  r
}

def <- run("default")
pos <- run("poscounts")

cat("DEFAULT  : tested", sum(!is.na(def$pvalue)),
    "| padj<0.05:", sum(def$padj < 0.05, na.rm = TRUE),
    "| stat>0 (preterm):", sum(def$stat > 0, na.rm = TRUE), "\n")
cat("POSCOUNTS: tested", sum(!is.na(pos$pvalue)),
    "| padj<0.05:", sum(pos$padj < 0.05, na.rm = TRUE), "\n")

# side by side, top 8 by the default run's padj
cmp <- tibble(
  feature  = def$feature,
  p_def    = def$pvalue,   padj_def = def$padj, stat_def = def$stat,
  p_pos    = pos$pvalue,   padj_pos = pos$padj, stat_pos = pos$stat
)
cat("\nTop 8 by the default run's padj:\n")
print(cmp |> arrange(padj_def) |> head(8) |>
      mutate(across(where(is.numeric), \(x) signif(x, 3))), width = 130)

cat("\nSpecies significant in the default run but not in poscounts:\n")
cat(sum(!is.na(cmp$padj_def) & cmp$padj_def < 0.05 &
        (is.na(cmp$padj_pos) | cmp$padj_pos >= 0.05)), "of",
    sum(!is.na(cmp$padj_def) & cmp$padj_def < 0.05), "\n")
