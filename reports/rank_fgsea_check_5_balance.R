# Notebook 05, Exercise 1, Step 4 - is the ranking balanced enough to trust fgsea?

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

genes_by_outcome <- read_csv(here("data", "genes_by_outcome.csv"), show_col_types = FALSE)

stat_vec <- genes_by_outcome$stat
names(stat_vec) <- genes_by_outcome$feature
stat_vec <- sort(stat_vec[!is.na(stat_vec)], decreasing = TRUE)

ranked <- abs(stat_vec)
pct_pos <- function(x) round(100 * mean(x > 0), 1)

cat("== Share of positive statistics (stat > 0 = higher in PRETERM) ==\n")

sig_ids <- genes_by_outcome$feature[genes_by_outcome$padj < 0.05]
sig <- stat_vec[names(stat_vec) %in% sig_ids]

q500 <- ranked >= quantile(ranked, 1 - 500 / length(ranked))
top500 <- stat_vec[q500]

summary_tab <- tibble(
  slice        = c("all genes", "significant (padj<0.05)", "top 500 by |stat|"),
  n            = c(length(stat_vec), length(sig), length(top500)),
  pct_positive = c(pct_pos(stat_vec), pct_pos(sig), pct_pos(top500))
)
print(summary_tab)

cat("\nmean stat overall:", round(mean(stat_vec), 3), "(0 means no overall tilt)\n")
