# plot_all567_circular.R ----------------------------------------------------------
#
# Circular phylogenetic tree of ALL 567 protein accessions from
# R/accession_description_567.csv (mixed genes: COI, allergens, proteases, ...).
#
# Tip labels are plain accession numbers. The six accessions listed below are
# drawn in their own colours; everything else is dark grey.
#
# NOTE: mixing different genes in one alignment is a similarity view, not a
# strict evolutionary tree - see the "Choices made" note when we deliver it.
#
# Run:  Rscript R/plot_all567_circular.R

suppressPackageStartupMessages({
  library(ape)
  if (requireNamespace("here", quietly = TRUE)) root <- here::here() else root <- "."
})

tree <- read.tree(file.path(root, "R", "all_567.nwk"))
tree$tip.label <- gsub("'", "", tree$tip.label)   # FastTree -quote adds these

# --- colours: the six accessions to highlight --------------------------------
highlight <- c(
  "AAS93675.1" = "red",
  "AAS93667.1" = "blue",
  "AAS93669.1" = "darkgreen",
  "AAO15606.1" = "purple",
  "AAO15613.1" = "orange",
  "CAL18242.1" = "darkcyan"
)

tip_col <- rep("gray35", length(tree$tip.label))
hit <- match(names(highlight), tree$tip.label)
tip_col[hit[!is.na(hit)]] <- highlight[which(!is.na(hit))]
cat("highlighted found in tree:", sum(!is.na(hit)), "of", length(highlight), "\n")
cat("tips:", Ntip(tree), "\n")

# --- plot ---------------------------------------------------------------------
# Draw as a CLADOGRAM (use.edge.length = FALSE): with mixed genes the branch
# lengths vary hugely (0 to ~3 substitutions/site) and type="fan" uses them as
# radii, which blew the layout apart. Equal spacing keeps the circle readable.
# The true branch lengths live in R/all_567.nwk if you ever need them.

# Two-pass draw with DIAGONAL (radial) labels.
# Pass 1: measure the tree's extent. Pass 2: redraw with widened limits.
# ape's fan mode writes labels horizontally, which smears 567 names together,
# so we draw the tree bare, fetch each tip's x/y from ape's hidden
# last_plot.phylo object, and rotate every label to lie along its radius.
measure_and_plot <- function(cex_lab = 0.14) {
  plot(tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
       edge.width = 0.5)
  u <- par("usr")
  ext <- max(abs(u[1:2]), abs(u[3:4]))
  r <- 1.10 * ext   # tips reach near the page edge: longer arc per label
  plot(tree, type = "fan", use.edge.length = FALSE, show.tip.label = FALSE,
       edge.width = 0.5, x.lim = c(-r, r), y.lim = c(-r, r))

  lp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  tx <- lp$xx[seq_len(Ntip(tree))]
  ty <- lp$yy[seq_len(Ntip(tree))]
  rad <- sqrt(tx^2 + ty^2)
  deg <- atan2(ty, tx) * 180 / pi
  flip <- deg > 90 | deg < -90          # left half: flip so text stays upright
  srt <- deg; srt[flip] <- deg[flip] + 180
  adj1 <- rep(0, Ntip(tree)); adj1[flip] <- 1
  off <- rad + 0.05 * ext               # nudge labels just past the tips
  px <- tx / rad * off; py <- ty / rad * off

  fontv <- rep(1, Ntip(tree))
  fontv[tree$tip.label %in% names(highlight)] <- 2   # bold the six
  cexv <- rep(cex_lab, Ntip(tree))
  cexv[tree$tip.label %in% names(highlight)] <- cex_lab * 1.6  # bigger too
  # text() accepts only ONE srt per call, so loop tip by tip
  for (i in seq_len(Ntip(tree))) {
    text(px[i], py[i], tree$tip.label[i], srt = srt[i],
         adj = c(adj1[i], 0.5), cex = cexv[i], col = tip_col[i],
         font = fontv[i], xpd = NA)
  }
}

png(file.path(root, "reports", "all567_tree_circular.png"),
    width = 3400, height = 3400, res = 400)
par(mar = c(0.5, 0.5, 1.5, 0.5))
measure_and_plot()
title(main = "Sarcoptes scabiei - similarity tree of all 567 accessions",
      cex.main = 0.9)
legend("bottomleft",
       legend = names(highlight),
       col = highlight, pch = 15, cex = 0.45, bty = "n", title = "Highlighted", y.intersp = 0.7)
dev.off()

pdf(file.path(root, "reports", "all567_tree_circular.pdf"), width = 8.5, height = 8.5)
par(mar = c(0.5, 0.5, 1.5, 0.5))
measure_and_plot()
title(main = "Sarcoptes scabiei - similarity tree of all 567 accessions",
      cex.main = 0.9)
legend("bottomleft", legend = names(highlight), col = highlight,
       pch = 15, cex = 0.45, bty = "n", title = "Highlighted", y.intersp = 0.7)
dev.off()

cat("Wrote reports/all567_tree_circular.png and .pdf\n")
