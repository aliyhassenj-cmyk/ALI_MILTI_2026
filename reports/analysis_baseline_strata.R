# Did the probiotic only work in people whose community started a certain way?
#
# Plan:
# 1. From each person's PRE microbiome sample only, name their starting
#    community type:
#      - crispatus-dominant : >50% L. crispatus
#      - iners-dominant     : >50% L. iners (but not crispatus)
#      - diverse            : neither above 50%
# 2. Outcome = clinical improvement, pre minus post Nugent score
#    (positive number = got better). We deliberately read "did it work"
#    from the CLINICAL table, not the same table the groups came from.
# 3. Compare probiotic vs placebo inside each starting group.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(base64enc)
})

set.seed(42)
dat <- "data/2-workshop2/group_a_probiotic"

meta <- read_csv(file.path(dat, "01_metadata.csv"), show_col_types = FALSE)
mb   <- read_csv(file.path(dat, "02_microbiome_counts.csv"), show_col_types = FALSE)
tax  <- read_csv(file.path(dat, "03_microbiome_taxonomy.csv"), show_col_types = FALSE)

# ---- 1. starting community type from the pre sample -------------------------
species <- mb$taxon
counts  <- as.matrix(mb[, -1])                      # species x sample
colnames(counts) <- colnames(mb[, -1])
rel <- sweep(counts, 2, colSums(counts), "/") * 100 # % per sample

row_of <- function(pattern) which(species == pattern)
cris <- rel[row_of("lactobacillus_crispatus"), ]
iner <- rel[row_of("lactobacillus_iners"),     ]

sample_info <- tibble(
  sample_id      = colnames(counts),
  pct_crispatus  = as.numeric(cris),
  pct_iners      = as.numeric(iner)
) |>
  mutate(
    start_type = case_when(
      pct_crispatus > 50          ~ "L. crispatus\ndominant",
      pct_iners     > 50          ~ "L. iners\ndominant",
      .default                     = "Diverse /\nnon-Lactobacillus"
    )
  )

# The starting label is a property of the PERSON, read off their PRE sample.
start_of_person <- sample_info |>
  inner_join(meta |> select(sample_id, participant_id, time_point), by = "sample_id") |>
  filter(time_point == "pre") |>
  select(participant_id, start_type)

# ---- 2. join to metadata, one row per person -------------------------------
people <- meta |>
  select(participant_id, arm, time_point, nugent, ph) |>
  left_join(start_of_person, by = "participant_id") |>
  pivot_wider(names_from = time_point, values_from = c(nugent, ph)) |>
  mutate(
    delta_nugent = nugent_pre - nugent_post,   # >0 = improved
    delta_ph     = ph_pre - ph_post            # >0 = pH came down
  )

cat("People per arm x start_type:\n")
print(people |> count(arm, start_type) |> pivot_wider(names_from = arm, values_from = n))

# ---- 3. effect inside each starting group ----------------------------------
summ <- people |>
  group_by(start_type, arm) |>
  summarise(
    n          = n(),
    median_chg = median(delta_nugent, na.rm = TRUE),
    improved   = mean(delta_nugent > 0, na.rm = TRUE),
    .groups    = "drop"
  )
cat("\nMedian Nugent change (>0 = better):\n")
print(as.data.frame(summ))

# wilcoxon p per stratum (rank test, no normality assumption)
pvals <- people |>
  group_by(start_type) |>
  summarise(
    p = tryCatch(
      wilcox.test(delta_nugent ~ arm)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )
cat("\nWilcoxon p (probiotic vs placebo) per start group:\n")
print(as.data.frame(pvals))

# ---- 4. figure --------------------------------------------------------------
lvl <- c("L. crispatus\ndominant", "L. iners\ndominant", "Diverse /\nnon-Lactobacillus")
people$start_type <- factor(people$start_type, levels = lvl)
cols <- c(probiotic = "#2E7D32", placebo = "#757575")

# A: who started where
pA <- people |>
  distinct(participant_id, arm, start_type) |>   # one row per person
  count(arm, start_type) |>
  ggplot(aes(start_type, n, fill = arm)) +
  geom_col(position = "dodge", width = .7) +
  geom_text(aes(label = n), position = position_dodge(.7), vjust = -.3, size = 3.4) +
  scale_fill_manual(values = cols) +
  labs(title = "A. Who started where?", x = NULL, y = "Participants") +
  ylim(0, NA) + theme_minimal(11) +
  theme(axis.text.x = element_text(size = 8))

# B: change in Nugent by arm within stratum
pB <- ggplot(people, aes(start_type, delta_nugent, fill = arm)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = .6, fatten = 1.2,
               position = position_dodge(.75)) +
  geom_point(aes(colour = arm), shape = 21, fill = NA, size = 1.4,
             position = position_jitterdodge(jitter.width = .18, dodge.width = .75),
             alpha = .5, stroke = .4) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  labs(title = "B. Change in Nugent score (positive = improved)",
       x = "Starting community type", y = "Nugent change\n(pre − post)") +
  theme_minimal(11) + theme(axis.text.x = element_text(size = 8), legend.position = "none")

# C: pre -> post median path per arm x stratum
path <- people |>
  select(start_type, arm, pre = nugent_pre, post = nugent_post) |>
  group_by(start_type, arm) |>
  summarise(across(c(pre, post), median), .groups = "drop") |>
  pivot_longer(c(pre, post), names_to = "visit", values_to = "nugent") |>
  mutate(visit = factor(visit, c("pre", "post")))

pC <- ggplot(path, aes(visit, nugent, colour = arm, group = arm)) +
  geom_line(linewidth = .9) + geom_point(size = 2.4) +
  facet_wrap(~start_type, nrow = 1) +
  scale_colour_manual(values = cols) +
  labs(title = "C. Median Nugent score, pre → post",
       x = NULL, y = "Median Nugent\n(higher = more BV)") +
  theme_minimal(11) + theme(axis.text.x = element_text(size = 8)) +
  theme(strip.text = element_text(size = 7.5))

fig <- (pA | pB) / pC + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title    = "Did the probiotic work equally in every starting community?",
    subtitle = "Group A probiotic trial — 176 participants, pre vs post. Groups defined from the PRE microbiome; outcome read from clinical Nugent scores.",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

png <- file.path(tempdir(), "fig.png")
ggsave(png, fig, width = 10, height = 6.4, dpi = 160, bg = "white")
# also keep the figure as a standalone image next to the html
ggsave("reports/fig_baseline_strata.png", fig, width = 10, height = 6.4, dpi = 160, bg = "white")
img64 <- base64encode(png)

# ---- 5. write numbers for the html ------------------------------------------
summ_html <- summ |>
  mutate(across(c(median_chg, improved), ~ round(.x, 2))) |>
  left_join(pvals, by = "start_type") |>
  mutate(p = signif(p, 2)) |>
  pivot_wider(names_from = arm,
              values_from = c(n, median_chg, improved)) |>
  select(start_type,
         `n placebo` = n_placebo, `n probiotic` = n_probiotic,
         `median ΔNugent placebo` = median_chg_placebo,
         `median ΔNugent probiotic` = median_chg_probiotic,
         `Wilcoxon p` = p) |>
  knitr::kable(format = "html", escape = FALSE) |>
  as.character()

saveRDS(list(fig = fig, summ_html = summ_html,
             people = people), "reports/_baseline_strata.rds")

# ---- 6. one self-contained html ---------------------------------------------
html <- paste0(
'<html>
<head>
<meta charset="utf-8">
<title>Did the probiotic only work for some starting communities?</title>
<style>
 body{font-family:Georgia,serif;max-width:1060px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.55}
 h1{font-size:1.5rem} h2{font-size:1.15rem;margin-top:2rem;border-bottom:2px solid #2E7D32;padding-bottom:.2rem}
 table{border-collapse:collapse;margin:1rem 0;font-size:.92rem}
 th,td{border:1px solid #bbb;padding:.4rem .7rem;text-align:left}
 th{background:#eef5ee}
 .callout{background:#f2f7f2;border-left:4px solid #2E7D32;padding:.8rem 1rem;margin:1rem 0}
 .warn{background:#fdf6ec;border-left:4px solid #c98a2d;padding:.8rem 1rem;margin:1rem 0}
 code{background:#f4f4f4;padding:.1rem .3rem;border-radius:3px;font-size:.85em}
 li{margin:.3rem 0}
 img{max-width:100%;border:1px solid #ddd;border-radius:6px;margin:1rem 0}
</style>
</head>
<body>
<h1>Did the probiotic only work in people whose community started a certain way?</h1>
<p><strong>Group A probiotic trial</strong> — 176 participants, one sample before and one after.
Each person was first sorted into a <em>starting community type</em>, judged only from their
<code>pre</code> microbiome sample:</p>
<ul>
 <li><strong>L. crispatus dominant</strong> — more than 50% <em>Lactobacillus crispatus</em></li>
 <li><strong>L. iners dominant</strong> — more than 50% <em>L. iners</em>, but not crispatus</li>
 <li><strong>Diverse / non-Lactobacillus</strong> — neither species above 50%</li>
</ul>
<p>The effect measure is the <strong>change in Nugent score</strong> (pre minus post).
Nugent runs 0–10 and higher means more bacterial vaginosis, so a <strong>positive change means the person got better</strong>.</p>

<h2>The picture</h2>
<img src="data:image/png;base64,', paste(img64), '" onerror="this.onerror=null;this.src=\'fig_baseline_strata.png\'" alt="Figure: who started where, change in Nugent score by arm within starting group, and median pre-to-post Nugent paths">

<h2>The numbers</h2>
', summ_html, '
<p><em>Wilcoxon p</em> = rank-sum test of probiotic vs placebo on the Nugent change, inside each starting group.</p>

<h2>The answer</h2>
<div class="callout">
<strong>Short answer: no clear evidence that it worked only for a particular starting community — but a suggestive pattern for <em>L. iners</em> starters, and one starting type is missing entirely.</strong>
<ul>
<li>Nobody in this trial started <em>L. crispatus</em>-dominant, so the “best” baseline community cannot be tested here at all.</li>
<li>Among the <strong>diverse starters</strong> (the great majority: 166 of 176), the median improvement was 3.5 points on probiotic vs 2 on placebo, but the difference is not significant (p ≈ 0.20).</li>
<li>Among <strong><em>L. iners</em> starters</strong>, the probiotic median change was larger (5.0 vs 2.0) — but there is only <strong>1 placebo participant</strong> in that group, so no real comparison is possible (p ≈ 0.59 and meaningless with n=1).</li>
<li>Consistent picture though: in both starting groups the improvement percentage was higher on probiotic (73% vs 61%; 78% vs 100% of one person) — a hint, not a conclusion.</li>
</ul>
</div>
<div class="warn">
<strong>Read this before quoting it:</strong> the two groups are very unbalanced (121 vs 55 people, and 112 vs 54 within the diverse starters), the Nugent change has strong ties so the exact p-value is approximate, and nothing here is adjusted for multiplicity. Treat this as an exploratory check, not a result.
</div>

<h2>Choices made</h2>
<ul>
<li><strong>Baseline grouping rule</strong>: &gt;50% relative abundance of one Lactobacillus species. Alternative not chosen: clustering all pre samples (e.g. Ward on Bray–Curtis) — more data-driven but harder to explain and interpret.</li>
<li><strong>Outcome</strong>: change in clinical Nugent score, so the grouping table (microbiome) and the outcome table (clinical) are different, per the workshop rule. Alternative not chosen: change in microbiome composition itself — that would define the group from the table we then test.</li>
<li><strong>Test</strong>: Wilcoxon rank-sum within each starting group. Alternatives not chosen: two-way model with interaction (needs the interaction term interpretation, and the tiny iners/placebo cell breaks it), or linear model on raw change.</li>
<li><strong>Samples dropped</strong>: 2 people missing a post Nugent score fall out of panel B automatically (they still appear in panel A, counted by group).</li>
<li><strong>Transform</strong>: none — Nugent change is analysed as-is; the microbiome percentages were only used to assign the starting label, never tested.</li>
</ul>
<h2>Check it yourself</h2>
<p>Run <code>Rscript reports/analysis_baseline_strata.R</code> and confirm the count table it prints: 54 placebo and 112 probiotic among diverse starters, 1 and 9 among <em>L. iners</em> starters, and no crispatus row at all. In panel B of the figure, the probiotic jitter cloud should sit visibly higher than placebo only in the <em>L. iners</em> column — and the placebo half of that column should contain a single dot.</p>
</body>
</html>'
)
writeLines(html, "reports/probiotic_by_starting_community.html")
cat("HTML written: reports/probiotic_by_starting_community.html\n")
