---
title: "Phylogenetic analysis of unique protein sequences from *Sarcoptes scabiei*"
---

# Phylogenetic analysis of unique protein sequences from *Sarcoptes scabiei*

## Figure 1

![Figure 1. Circular phylogenetic tree of the 93 unique protein sequences from *Sarcoptes scabiei*.](reports/all_unique_tree_circular.png)

**Figure 1.** Circular phylogenetic tree of the 93 unique protein sequences retained after
removal of duplicates among the 567 retrieved accessions (each tip label shows the
representative accession and, in parentheses, the number of accessions it represents).
The tree was inferred with FastTree 2.1.11 from an amino-acid alignment using the
Jones–Taylor–Thornton (JTT) model with a CAT approximation of 20 rate categories;
branch support was assessed with the Shimodaira–Hasegawa-like test (1,000 resamples)
and is stored in the tree file. Tips are drawn at equal branch lengths for readability;
the underlying branch lengths, expressed as mean expected amino-acid substitutions per
site (range 0–1.46; total tree length 14.7), are provided in the Newick file
(`all_unique.nwk`). The six accessions of interest for this study are highlighted in
colour; all other tips are shown in dark grey. The number in parentheses after each
accession, (n=k), gives the number of the 567 retrieved accessions that encode an
identical sequence and are thus represented by that tip.

## Colour scheme of Figure 1

Two kinds of colour are used in Figure 1, each with a different purpose:

1. **Dark grey for the background tree (87 of 93 tips).** Grey is used deliberately so
   that the surrounding sequences recede visually and act as phylogenetic context
   rather than the message of the figure. A neutral grey also prints correctly in both
   colour and greyscale.

2. **Six high-contrast colours for the accessions of interest.** The six accessions
   that this study focuses on were each assigned one distinct colour — red
   (AAS93675.1), blue (AAS93667.1), dark green (AAS93669.1), purple (AAO15606.1),
   orange (AAO15613.1) and dark cyan (CAL18242.1). These colours were selected
   because (i) they are maximally distinguishable from one another, so each focus
   sequence can be traced individually around the circular layout, and (ii) they
   differ strongly in both hue and lightness from the grey background, so the focus
   tips stand out immediately even at the small label size required by a 93-tip
   circular tree. Bold font was applied to the coloured labels for the same reason.

In manuscript form, this colour scheme lets the reader answer the study's key
question — where the six focus sequences sit relative to the rest of the *S. scabiei*
proteins — at a glance, without having to search 93 labels.

## Interpretation of the tip labels (n = …)

Each tip label in Figure 1 has the form **ACCESSION (n=k)**, where the two parts have
the following meaning:

- **ACCESSION** is the GenBank accession of the representative sequence, i.e. the
  one sequence from each group of identical sequences that was actually included in
  the alignment and the tree.
- **(n=k)** gives the number of accessions that this representative stands for. Of
  the 567 accessions originally retrieved from NCBI, many proved to be exactly the
  same protein sequence (identical submissions from independent studies or repeated
  uploads). All identical sequences were therefore collapsed into a single unique
  sequence, and *k* is the size of that group.

For example, the tip **WHU50584.1 (n=34)** indicates that 34 of the 567 retrieved
accessions encode this identical cytochrome c oxidase subunit I sequence, all of
which are represented by this single tip. A tip labelled **(n=1)** denotes an
accession whose sequence was unique among the retrieved set — as is the case for
most of the non-COI proteins (myosin, allergens, metabolic enzymes, etc.).

Collapsing identical sequences before tree inference is standard practice: strictly
identical sequences contribute no phylogenetic information (they are joined by
branches of length zero) but would add 474 unreadable labels to the figure. The
*n* values over all 93 tips sum to 567, so every retrieved accession is accounted
for in the tree. Because duplicate sequences carry no information about the tree
shape, the counts in parentheses are labels for readability only and do not affect
branch lengths or topology.

## Table 1. Description and source of each unique protein

All 93 unique protein sequences were identified against the NCBI protein database
(GenBank). Every sequence originates from *Sarcoptes scabiei*. The full list is given
below; representative accessions are those shown as tip labels in Figure 1.


| Representative accession | Protein description (NCBI) | Source organism |
|---|---|---|
| BAJ10006.1 | myosin heavy chain partial | *Sarcoptes scabiei* |
| XUS00055.1 | trypsin-like serine protease partial | *Sarcoptes scabiei* |
| WHU50584.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50563.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50542.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50537.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50505.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50500.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WHU50462.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WOP59345.1 | voltage-sensitive sodium channel partial | *Sarcoptes scabiei* |
| WKB12265.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| WKB12215.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| UXI21683.1 | enoyl-CoA hydratase domain-containing protein 3 | *Sarcoptes scabiei* |
| UXI20073.1 | isobutyryl-CoA dehydrogenase | *Sarcoptes scabiei* |
| UXI19041.1 | small nuclear ribonucleoprotein-associated protein B-like | *Sarcoptes scabiei* |
| UXI17943.1 | methylcrotonoyl-CoA carboxylase beta chain | *Sarcoptes scabiei* |
| UXI16844.1 | heat shock protein HSP 90-alpha | *Sarcoptes scabiei* |
| KAI7693147.1 | Vacuolar protein sorting-associated protein 29 | *Sarcoptes scabiei* |
| KAI7691166.1 | Tetratricopeptide repeat protein 4 | *Sarcoptes scabiei* |
| KAI7690134.1 | hypothetical protein SSS 08815 | *Sarcoptes scabiei* |
| KAI7690133.1 | hypothetical protein SSS 08815 | *Sarcoptes scabiei* |
| KAI7690132.1 | Glutathione-specific gamma-glutamylcyclotransferase 1 | *Sarcoptes scabiei* |
| KAF7495932.1 | Putative cytosolic Fe-S cluster assembly factor | *Sarcoptes scabiei* |
| KAF7494399.1 | Double-strand-break repair protein rad21 -like protein | *Sarcoptes scabiei* |
| KAF7491894.1 | Triosephosphate isomerase | *Sarcoptes scabiei* |
| KAF7491437.1 | Xylulose kinase | *Sarcoptes scabiei* |
| KAF7487801.1 | Glutathione-specific gamma-glutamylcyclotransferase 1 | *Sarcoptes scabiei* |
| UHQ62601.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| UBP77630.1 | cytochrome c oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02212.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02210.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02195.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02194.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02192.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02189.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02183.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02161.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02153.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QIK02152.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| LAB77524.1 | xylulose kinase-like protein | *Sarcoptes scabiei* |
| LAB77528.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB85655.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB85663.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB98550.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB74765.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB80639.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB83486.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB83502.1 | hypothetical protein | *Sarcoptes scabiei* |
| LAB97190.1 | Calcineurin-like phosphoesterase domain containing protein partial partial | *Sarcoptes scabiei* |
| QNG98905.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70654.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70644.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70628.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70612.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70609.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70564.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70554.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70526.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| QCC70524.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AVW83288.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AVW83287.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AVW83275.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AAO15613.1 | major allergen 1 partial | *Sarcoptes scabiei* |
| AAO15612.1 | paramyosin partial | *Sarcoptes scabiei* |
| AAO15611.1 | vitellogenin-like protein partial | *Sarcoptes scabiei* |
| AAO15610.1 | cytochrome Bc1 complex chain B-like protein partial | *Sarcoptes scabiei* |
| AAO15609.1 | cytochrome oxidase subunit II partial (mitochondrion) | *Sarcoptes scabiei* |
| AAO15608.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AAO15606.1 | cathepsin L-like protease partial | *Sarcoptes scabiei* |
| KPM10469.1 | Phosphorylated CTD-interacting factor-like protein | *Sarcoptes scabiei* |
| KPM08840.1 | Calcineurin-like phosphoesterase domain containing protein partial | *Sarcoptes scabiei* |
| KPM05597.1 | double-strand-break repair protein rad21-like protein | *Sarcoptes scabiei* |
| KPM02699.1 | tetratricopeptide repeat containing cytosolic Fe-S cluster assembly factor-like protein | *Sarcoptes scabiei* |
| ALA99922.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| ALA99921.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| ALA99917.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AKG58795.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AIB63974.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AIB63972.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AIB63970.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| AIB63969.1 | cytochrome oxidase subunit I partial (mitochondrion) | *Sarcoptes scabiei* |
| CAL18242.1 | hypothetical protein partial | *Sarcoptes scabiei* |
| AAO15607.1 | glutathione S-transferase | *Sarcoptes scabiei* |
| AAS93672.1 | Sar s 1 allergen SMIPP-C Yv4025A02 | *Sarcoptes scabiei* |
| AAS93673.1 | Sar s 1 allergen SMIPP-C Yv6008G08 | *Sarcoptes scabiei* |
| AAS93675.1 | Sar s 1 allergen SMIPP-C Yv5009F04 | *Sarcoptes scabiei* |
| AAS93676.1 | Sar s 1 allergen SMIPP-C Yv6018B11 partial | *Sarcoptes scabiei* |
| AAS93674.1 | Sar s 1 allergen SMIPP-C Yv4028C12 | *Sarcoptes scabiei* |
| AAS93667.1 | Sar s 1 allergen Yv4003H01 | *Sarcoptes scabiei* |
| AAS93668.1 | Sar s 1 allergen Yv5020C01 | *Sarcoptes scabiei* |
| AAS93670.1 | Sar s 1 allergen Yv6030H07 | *Sarcoptes scabiei* |
| AAS93669.1 | Sar s 1 allergen Yv5032C08 | *Sarcoptes scabiei* |
| AAS93671.1 | Sar s 1 allergen Yv9053H09 partial | *Sarcoptes scabiei* |


**Table 1 (summary).** The 93 unique proteins fall into the following groups:

| Protein group | Number of unique sequences |
|---|---|

| Cytochrome c oxidase subunit I (COI, mitochondrial) | 44 |
| Sar s 1 allergens / SMIPP-C variants | 11 |
| Hypothetical proteins (unknown function) | 11 |
| Structural muscle proteins (myosin, paramyosin) | 2 |
| Proteases (trypsin-like, cathepsin L-like) | 2 |
| Metabolic enzymes | 6 |
| Other (GST, HSP90, vitellogenin-like, sodium channel, rad21-like, etc.) | 17 |
| **Total** | **93** |

## Methods (tree construction and scale)

Protein sequences were aligned and a phylogeny was inferred with FastTree 2.1.11
(double precision) using initial BLOSUM45 distances, the JTT substitution model with
a CAT approximation of 20 rate categories, and SH-like support values from 1,000
resamples. Branch lengths are in mean expected amino-acid substitutions per site;
a scale bar of **0.1 substitutions per site** is recommended for any redrawn version
of the tree that uses real branch lengths. Support values range from 0 to 1 and are
stored at each internal node in `all_unique.nwk`.

## Choices made

- Tree shown as an equal-length (cladogram) circular layout; alternative not chosen:
  plotting real branch lengths, which are too uneven (0–1.46) to use as radii in a
  circular figure.
- One tip per unique sequence (93 tips) with duplicate counts in the label;
  alternative not chosen: all 567 accessions, which is unreadable.
- Background tips in dark grey; alternative not chosen: black, which competes with
  the coloured focus tips.
- Recommended scale bar of 0.1 amino-acid substitutions per site; alternative not
  chosen: 0.05 or 0.5, which would be too small or too large relative to the
  observed branch lengths.
- Note: the published circular figure does not display a scale bar because branch
  lengths were not drawn; the scale statement belongs in the caption (as above) or
  in a redrawn version.
- Identical sequences (567 accessions) were collapsed to 93 unique representatives
  for tree inference, with duplicate counts shown as (n=k) on each tip;
  alternative not chosen: inferring the tree from all 567 accessions, which wastes
  computation and produces an unreadable figure without changing the topology.
