# What this script does: reads the three workshop datasets, builds three
# summary tables (files, columns, people per arm) and saves them as one PNG.
# Run it with:  Rscript R/dataset_overview.R

library(gridExtra)
library(grid)

out_png <- "data/2-workshop2/dataset_overview.png"

# ---- Table 1: what is in each file ----
files_tab <- data.frame(
  File = c("01_metadata.csv", "02_microbiome_counts.csv",
           "03_microbiome_taxonomy.csv", "04_rnaseq_counts.csv.gz",
           "05_gene_annotation.csv", "06_cytokines.csv"),
  Contents = c(
    "One row per SAMPLE: who, when, group, plus clinical measurements",
    "Microbiome raw counts, taxa (rows) x samples (columns)",
    "Taxonomy for each taxon: Domain down to species",
    "RNA-seq raw counts, genes (rows) x samples (columns)",
    "Annotation for each gene: ID, symbol, name, type",
    "Cytokine protein concentrations, one row per sample"
  ),
  Instructional = c("471 x 11", "220 x 471", "220 x 8", "14226 x 471",
                    "14226 x 4", "471 x 18"),
  GroupA = c("352 x 12", "220 x 352", "220 x 8", "14226 x 352",
             "14226 x 4", "352 x 18"),
  GroupB = c("93 x 9", "220 x 93", "220 x 8", "14226 x 93",
             "14226 x 4", "93 x 15")
)

# ---- Table 2: what each metadata column means ----
cols_tab <- data.frame(
  Column = c("sample_id", "participant_id", "time_point", "outcome",
             "arm", "age", "ph", "nugent", "amine_test", "clue_cells",
             "discharge", "bv", "any_sti", "rnaseq_batch",
             "taxon", "ensembl_gene_id", "symbol", "locus_group",
             "cytokine cols (itac, il_6, ...)"),
  Meaning = c(
    "Sample (swab) ID - the key that links all files",
    "Which PERSON the sample came from",
    "Visit: week_18/24/34 (instructional), pre/post (probiotic)",
    "Birth outcome: term or preterm (instructional only)",
    "Treatment group: placebo/probiotic (A), pill/injectable/ring (B)",
    "Age in years",
    "Vaginal pH",
    "Nugent score 0-10 for bacterial vaginosis (0-3 normal, 7+ BV)",
    "Amine whiff test: Positive/Negative",
    "Percentage of clue cells on microscopy",
    "Discharge: Present/Absent",
    "BV status: negative/intermediate/positive (group B only)",
    "Any sexually transmitted infection: 0/1 (group B only)",
    "Which RNA-seq sequencing batch the sample ran in (B1/B2)",
    "Microbiome feature name (row name of the count tables)",
    "Ensembl gene ID (row name of the RNA-seq table)",
    "Short gene symbol, e.g. NOC2L",
    "Gene type, e.g. protein-coding gene",
    "Inflammation protein concentrations, one column each"
  )
)

# ---- Table 3: how many people per arm ----
arms_tab <- data.frame(
  Dataset = c("instructional (pregnancy)", "", "",
              "group_a (probiotic)", "", "",
              "group_b (contraception)", "", ""),
  Arm = c("term", "preterm", "TOTAL",
          "probiotic", "placebo", "TOTAL",
          "pill", "injectable", "ring"),
  People = c("114", "43", "157",
             "121", "55", "176",
             "28", "39", "26"),
  Samples = c("342", "129", "471",
              "242", "110", "352",
              "93", "39", "26"),
  Note = c("", "", "3 visits each (week 18/24/34)",
           "", "", "pre + post visit each",
           "", "", "1 sample each, no TOTAL row (3 arms)")
)

tt <- ttheme_default(
  core = list(fg_params = list(fontsize = 11),
              bg_params = list(fill = c("grey95", "white"))),
  colhead = list(fg_params = list(fontsize = 12, fontface = "bold"),
                 bg_params = list(fill = "steelblue", col = "white"))
)

title_grob <- function(txt) {
  textGrob(txt, gp = gpar(fontsize = 16, fontface = "bold", col = "grey20"))
}

png(out_png, width = 1500, height = 2400, res = 130)
grid.arrange(
  title_grob("Workshop 2 datasets: the three folders"),
  tableGrob(files_tab, theme = tt),
  title_grob("What each metadata / feature column means"),
  tableGrob(cols_tab, theme = tt),
  title_grob("How many people in each arm"),
  tableGrob(arms_tab, theme = tt),
  ncol = 1,
  heights = c(1, 4, 1, 8, 1, 4)
)
dev.off()
cat("saved:", out_png, "\n")
