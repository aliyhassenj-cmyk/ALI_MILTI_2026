# plot_coi_tree.R ----------------------------------------------------------------
#
# Draws the phylogenetic tree of the 43 unique COI (cytochrome oxidase
# subunit I) protein sequences from R/coi_accession_map.csv + R/coi_unique.nwk.
#
# The raw table had 508 accessions, but FastTree showed only 43 unique
# sequences. Each tip is labelled "REPRESENTATIVE ACCESSION (n=HOW MANY)".
#
# Run:  Rscript R/plot_coi_tree.R

suppressPackageStartupMessages({
  library(readr)
  library(ape)
  if (requireNamespace("here", quietly = TRUE)) root <- here::here() else root <- "."
})

tree <- read.tree(file.path(root, "R", "coi_unique.nwk"))
map  <- read_csv(file.path(root, "R", "coi_accession_map.csv"), show_col_types = FALSE)

# FastTree -quote wraps names in single quotes; strip them or nothing matches
tree$tip.label <- gsub("'", "", tree$tip.label)

n_per_rep <- table(map$rep)
tree$tip.label <- paste0(tree$tip.label, "  (n=", n_per_rep[tree$tip.label], ")")

cat("tips:", Ntip(tree), "\n")

png(file.path(root, "reports", "coi_tree.png"), width = 2200, height = 2000, res = 200)
par(mar = c(4, 1, 3, 1))
plot(tree, type = "phylogram", show.tip.label = TRUE, cex = 0.7,
     label.offset = 0.003, font = 1, edge.width = 1.2)
axisPhylo(1, las = 1)
title(main = "Sarcoptes scabiei COI - FastTree of 43 unique sequences",
      cex.main = 1.0, line = 1)
mtext("branch length = substitutions per site (protein)", side = 1, line = 2.2, cex = 0.75)
dev.off()

pdf(file.path(root, "reports", "coi_tree.pdf"), width = 10, height = 9)
par(mar = c(4, 1, 3, 1))
plot(tree, type = "phylogram", show.tip.label = TRUE, cex = 0.55,
     label.offset = 0.003, font = 1, edge.width = 1.2)
axisPhylo(1, las = 1)
title(main = "Sarcoptes scabiei COI - FastTree of 43 unique sequences")
mtext("branch length = substitutions per site (protein)", side = 1, line = 2.2, cex = 0.75)
dev.off()

cat("Wrote reports/coi_tree.png and reports/coi_tree.pdf\n")
