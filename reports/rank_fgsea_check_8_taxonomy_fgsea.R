# Notebook 05, Exercise 2, Steps 1-3 - taxonomic sets and fgsea on species

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(fgsea)
})

species_by_outcome <- read_csv(here("data", "species_by_outcome.csv"), show_col_types = FALSE)

# ---- the ranking: stat named by species, NAs dropped --------------------------
species_stat <- species_by_outcome$stat
names(species_stat) <- species_by_outcome$feature
species_stat <- sort(species_stat[!is.na(species_stat)], decreasing = TRUE)
cat("species in ranking:", length(species_stat), "\n")

# ---- identifier checks --------------------------------------------------------
cat("species missing a Genus:",
    sum(is.na(species_by_outcome$Genus[match(names(species_stat), species_by_outcome$feature)])), "\n")
cat("species in two sets:",
    sum(duplicated(species_by_outcome$feature)), "(0 means one genus per species)\n")

# ---- Step 1: one set per genus, at least two species --------------------------
genera <- species_by_outcome |> distinct(feature, Genus)
genus_counts <- table(genera$Genus)
genus_sets <- genera |>
  filter(Genus %in% names(genus_counts[genus_counts >= 2])) |>
  with(split(feature, Genus))
cat("genus sets:", length(genus_sets),
    "| sizes:", paste(sort(lengths(genus_sets), decreasing = TRUE), collapse = ", "), "\n")

fg <- fgsea(stats = species_stat, pathways = genus_sets, minSize = 2, maxSize = 50)
fg <- as.data.frame(fg) |> arrange(padj) |> as_tibble()
cat("\n== Step 1: fgsea on genus sets ==\n")
print(fg |> select(pathway, size, NES, pval, padj) |>
      mutate(across(where(is.numeric), \(x) signif(x, 3))), n = 10)
cat("any padj < 0.05:", any(fg$padj < 0.05), "\n")

# ---- Step 2: family instead of genus ------------------------------------------
family_counts <- table(species_by_outcome$Family)
family_sets <- species_by_outcome |>
  filter(Family %in% names(family_counts[family_counts >= 2])) |>
  distinct(feature, Family) |>
  with(split(feature, Family))
cat("\n== Step 2: family sets ==")
cat("| sets:", length(family_sets), "| sizes:", paste(lengths(family_sets), collapse = ", "), "\n")
fg_fam <- fgsea(stats = species_stat, pathways = family_sets, minSize = 2, maxSize = 50)
fg_fam <- as.data.frame(fg_fam) |> arrange(padj) |> as_tibble()
print(fg_fam |> select(pathway, size, NES, pval, padj) |>
      mutate(across(where(is.numeric), \(x) signif(x, 3))), n = 10)

top_genus  <- fg$pathway[1]
top_family <- fg_fam$pathway[1]
cat("\ntop genus:", top_genus, "- its family:",
    species_by_outcome$Family[species_by_outcome$Genus == top_genus][1], "\n")
cat("top family:", top_family, "\n")
cat("top-genus species inside the top family:",
    paste(intersect(genus_sets[[top_genus]], family_sets[[top_family]]), collapse = ", "), "\n")

# ---- Step 3: the simple way - the 53 species sorted by stat -------------------
out <- species_by_outcome |>
  select(feature, Genus, stat, padj) |>
  arrange(desc(stat))
write_csv(out, here("data", "species_sorted_by_stat.csv"))
cat("\n== Step 3: species sorted by stat (top 12) ==\n")
print(out |> head(12) |> mutate(across(where(is.numeric), \(x) signif(x, 3))), n = 12)
