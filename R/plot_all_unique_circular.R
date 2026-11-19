# plot_all_unique_circular.R -------------------------------------------------------
#
# SMALL, READABLE version of the diagonal accession tree.
#
# The 567 accessions collapse to 93 UNIQUE protein sequences, so instead of
# 567 unreadable labels we draw one tip per unique sequence. Each tip label is
# "REPRESENTATIVE ACCESSION (n=HOW MANY DUPLICATES)".
#
# Inputs:  R/all_unique.nwk          FastTree of the 93 unique sequences
#          R/all_accession_map.csv   accession -> representative mapping
# Output:  reports/all_unique_tree_circular.png / .pdf
#
# Run:  Rscript R/plot_all_unique_circular.R

suppressPackageStartupMessages({
  library(readr)
  library(ape)
  if (requireNamespace("here", quietly = TRUE)) root <- here::here() else root <- "."
})

tree <- read.tree(file.path(root, "R", "all_unique.nwk"))
map  <- read_csv(file.path(root, "R", "all_accession_map.csv"), show_col_types = FALSE)
tree$tip.label <- gsub("'", "", tree$tip.label)   # FastTree -quote adds these

# --- colours: the six accessions of interest ----------------------------------
highlight <- c(
  "AAS93675.1" = "red",
  "AAS93667.1" = "blue",
  "AAS93669.1" = "darkgreen",
  "AAO15606.1" = "purple",
  "AAO15613.1" = "orange",
  "CAL18242.1" = "darkcyan"
)

rep_of <- setNames(map$rep, map$accession)
hl_reps <- unname(rep_of[names(highlight)])
hit <- match(hl_reps, tree$tip.label)
cat("highlighted reps found in tree:", sum(!is.na(hit)), "of", length(hit), "\n")
cat("tips:", Ntip(tree), "\n")

# --- label each tip with how many accessions it represents --------------------
n_per_rep <- table(map$rep)
tree$tip.label <- paste0(tree$tip.label, " (n=", n_per_rep[tree$tip.label], ")")

tip_col <- rep("gray25", length(tree$tip.label))
tip_col[hit[!is.na(hit)]] <- highlight[!is.na(hit)]

# --- draw: circular tree, labels rotated to lie along each radius -------------
# use.edge.length = FALSE: mixing genes with very different lengths makes the
# real branch lengths unusable as radii; equal spacing keeps the circle readable.
draw <- function(cex_lab = 0.42) {
  plot(tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
       edge.width = 0.8)
  u <- par("usr"); ext <- max(abs(u[1:2]), abs(u[3:4]))
  plot(tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
       edge.width = 0.8, x.lim = c(-1.25 * ext, 1.25 * ext),
       y.lim = c(-1.25 * ext, 1.25 * ext))

  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  tx <- lp$xx[seq_len(Ntip(tree))]; ty <- lp$yy[seq_len(Ntip(tree))]
  rad <- sqrt(tx^2 + ty^2)
  deg <- atan2(ty, tx) * 180 / pi
  flip <- deg > 90 | deg < -90          # left half: flip so text stays upright
  srt <- deg; srt[flip] <- deg[flip] + 180
  adj1 <- rep(0, Ntip(tree)); adj1[flip] <- 1
  off <- rad + 0.06 * ext
  px <- tx / rad * off; py <- ty / rad * off

  fontv <- rep(1, Ntip(tree)); fontv[!is.na(hit)] <- 2          # bold the six
  cexv  <- rep(cex_lab, Ntip(tree)); cexv[!is.na(hit)] <- cex_lab * 1.4
  for (i in seq_len(Ntip(tree))) {
    text(px[i], py[i], tree$tip.label[i], srt = srt[i],
         adj = c(adj1[i], 0.5), cex = cexv[i], col = tip_col[i],
         font = fontv[i], xpd = NA)
  }
}

png(file.path(root, "reports", "all_unique_tree_circular.png"),
    width = 3000, height = 3000, res = 300)
par(mar = c(0.5, 0.5, 1.5, 0.5))
draw()
title(main = "Sarcoptes scabiei - 93 unique sequences (from 567 accessions)",
      cex.main = 1.0)
legend("bottomleft", legend = paste0(names(highlight), " -> ", hl_reps),
       col = highlight, pch = 15, cex = 0.8, bty = "n", title = "Highlighted",
       y.intersp = 0.8)
dev.off()

pdf(file.path(root, "reports", "all_unique_tree_circular.pdf"),
    width = 10, height = 10)
par(mar = c(0.5, 0.5, 1.5, 0.5))
draw()
title(main = "Sarcoptes scabiei - 93 unique sequences (from 567 accessions)",
      cex.main = 1.0)
legend("bottomleft", legend = paste0(names(highlight), " -> ", hl_reps),
       col = highlight, pch = 15, cex = 0.8, bty = "n", title = "Highlighted",
       y.intersp = 0.8)
dev.off()

cat("Wrote reports/all_unique_tree_circular.png and .pdf\n")
