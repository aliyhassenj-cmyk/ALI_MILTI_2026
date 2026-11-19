# plot_taxonomy_tree.R ----------------------------------------------------------
#
# Builds a taxonomy-based phylogenetic tree from 03_microbiome_taxonomy.csv.
#
# "Taxonomy-based" means the tree shape comes only from the ranks written in
# the table (Domain > Phylum > ... > species). Two species share a branch up
# to the rank where their labels stop agreeing. No DNA sequences are used, so
# branch LENGTHS here are just for drawing, not evolution time.
#
# Output (in reports/):
#   phylo_tree_genus.png / .pdf  - one tip per genus, tips coloured by phylum
#   phylo_tree_species.pdf       - one tip per species, taller page
#
# Run:  Rscript R/plot_taxonomy_tree.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ape)
  if (requireNamespace("here", quietly = TRUE)) {
    root <- here::here()
  } else {
    root <- "."
  }
})

dataset <- "instructional"
tax_file <- file.path(root, "data/2-workshop2/instructional_dataset",
                      "03_microbiome_taxonomy.csv")
tax <- read_csv(tax_file, show_col_types = FALSE, progress = FALSE)

# Fill in missing ranks so every lineage reaches the species label
# (microViz::tax_fix does the same job; we avoid the dependency here).
ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "species")
tax <- tax[ranks]
for (i in seq_along(ranks)[-1]) {
  miss <- is.na(tax[[i]]) | tax[[i]] == ""
  tax[[i]][miss] <- paste0(tax[[i - 1]][miss], "_unclassified")
}

# ---- lineage text -> Newick string -------------------------------------------
# A tiny nested list works like a folder tree: insert each lineage as a path,
# then walk the list to print Newick.
# R copies lists on assignment, so this function must RETURN the updated tree
# at every level - editing a copy in place would silently lose everything.
insert_path <- function(tree, path) {
  label <- path[1]
  child <- if (is.null(tree[[label]])) list() else tree[[label]]
  if (length(path) > 1) child <- insert_path(child, path[-1])
  tree[[label]] <- child
  tree
}

to_newick <- function(tree) {
  if (!length(tree)) return("")
  parts <- character(0)
  for (nm in names(tree)) {
    child <- to_newick(tree[[nm]])
    if (nzchar(child)) {
      parts <- c(parts, paste0(nm, "(", child, ")"))
    } else {
      parts <- c(parts, nm)  # a tip
    }
  }
  paste(parts, collapse = ",")
}

build_tree <- function(df, tip_rank) {
  keep <- ranks[1:which(ranks == tip_rank)]
  # Prefix every label with its rank (p_ = Phylum, g_ = Genus, ...). Without
  # this, labels that repeat across ranks - e.g. Saccharibacteria is BOTH a
  # phylum and a genus in this table - would collide into one node and the
  # tree would quietly gain duplicate tips.
  pref <- substr(keep, 1, 1)
  tree <- list()
  for (i in seq_len(nrow(df))) {
    path <- unique(paste0(pref, "_", vapply(keep, function(r) df[[r]][i], character(1))))
    tree <- insert_path(tree, path)
  }
  newick <- paste0("(", to_newick(tree), ");")
  phy <- read.tree(text = newick)
  phy$tip.label <- sub("^[a-z]_", "", phy$tip.label)  # back to clean names
  phy$edge.length <- rep(1, nrow(phy$edge))  # equal lengths: a cladogram
  phy
}

# ---- genus-level tree (the readable one) ------------------------------------
phy_g <- build_tree(tax, tip_rank = "Genus")

# phylum of each genus tip = the phylum of any species in that genus
genus_phylum <- tax |>
  distinct(Genus, Phylum) |>
  deframe()

n_phyla <- length(unique(genus_phylum[phy_g$tip.label]))
pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
         "#A65628", "#F781BF", "#66C2A5", "#8DA0CB", "#E5C494")
tip_cols <- pal[match(genus_phylum[phy_g$tip.label], unique(genus_phylum[phy_g$tip.label]))]

dir.create(file.path(root, "reports"), showWarnings = FALSE)

png(file.path(root, "reports", "phylo_tree_genus.png"),
    width = 2600, height = 4200, res = 220)
par(mar = c(1, 1, 3, 1))
plot(phy_g, type = "phylogram", use.edge.length = FALSE,
     show.tip.label = TRUE, cex = 0.65, font = 1,
     tip.color = tip_cols, label.offset = 0.6, no.margin = FALSE)
title(main = paste0("Taxonomy-based tree - ", dataset,
                    " (", phy_g$Nnode + 1, " genera, ", n_phyla, " phyla)"),
      cex.main = 1.1)
# a crude legend: one entry per phylum actually present
u_ph <- unique(genus_phylum[phy_g$tip.label])
legend("bottomright", legend = u_ph, col = pal[seq_along(u_ph)],
       pch = 15, cex = 0.8, bty = "n", title = "Phylum")
dev.off()

pdf(file.path(root, "reports", "phylo_tree_genus.pdf"), width = 12, height = 18)
par(mar = c(1, 1, 3, 1))
plot(phy_g, type = "phylogram", use.edge.length = FALSE,
     show.tip.label = TRUE, cex = 0.65, font = 1,
     tip.color = tip_cols, label.offset = 0.6)
title(main = paste0("Taxonomy-based tree - ", dataset, " (genera)"))
dev.off()

# ---- species-level tree (the tall one, saved as PDF only) --------------------
phy_s <- build_tree(tax, tip_rank = "species")
tip_cols_s <- pal[match(tax$Phylum[match(phy_s$tip.label, tax$species)],
                        unique(tax$Phylum))]
pdf(file.path(root, "reports", "phylo_tree_species.pdf"), width = 12, height = 30)
par(mar = c(1, 1, 3, 1))
plot(phy_s, type = "phylogram", use.edge.length = FALSE,
     show.tip.label = TRUE, cex = 0.55, font = 1,
     tip.color = tip_cols_s, label.offset = 0.6)
title(main = paste0("Taxonomy-based tree - ", dataset,
                    " (", phy_s$Nnode + 1, " species)"))
dev.off()

cat("Done.\n")
cat("Genera tips :", phy_g$Nnode + 1, "| unique genera in table:",
    length(unique(tax$Genus)), "\n")
cat("Species tips:", phy_s$Nnode + 1, "| species rows in table:", nrow(tax), "\n")
cat("Wrote reports/phylo_tree_genus.png, phylo_tree_genus.pdf, phylo_tree_species.pdf\n")
