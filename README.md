# 🧀 mi_cheese
Repository with code & data for the publication *Microbial interactions shape cheese flavour formation*. 

## 🖋️ Citation

If you use this resource, please cite:

 > Melkonian, C., Zorrilla, F., Kjærbølling, I. et al. Microbial interactions shape cheese flavour formation. Nat Commun 14, 8348 (2023). https://doi.org/10.1038/s41467-023-41059-2

## 👨🏻‍🍳 Abstract
Cheese fermentation and flavour formation are the result of complex biochemical reactions driven by the activity of multiple microorganisms. Here, we studied the roles of microbial interactions in flavour formation in a year-long Cheddar cheese making process, using a commercial starter culture containing Streptococcus thermophilus and Lactococcus strains. By using an experimental strategy whereby certain strains were left out from the starter culture, we show that S. thermophilus has a crucial role in boosting Lactococcus growth and shaping flavour compound profile. Controlled milk fermentations with systematic exclusion of single Lactococcus strains, combined with genomics, genome-scale metabolic modelling, and metatranscriptomics, indicated that S. thermophilus proteolytic activity relieves nitrogen limitation for Lactococcus and boosts de novo nucleotide biosynthesis. While S. thermophilus had large contribution to the flavour profile, Lactococcus cremoris also played a role by limiting diacetyl and acetoin formation, which otherwise results in an off-flavour when in excess. This off-flavour control could be attributed to the metabolic re-routing of citrate by L. cremoris from diacetyl and acetoin towards α-ketoglutarate. Further, closely related Lactococcus lactis strains exhibited different interaction patterns with S. thermophilus, highlighting the significance of strain specificity in cheese making. Our results highlight the crucial roles of competitive and cooperative microbial interactions in shaping cheese flavour profile.

## 🪤 Usage & description
Here is availalbe the majority of the processed data, including metabolomics, metatranscriptomics count tables and the corresponding code to reproduce the figures/results of the manuscript. The SLAB genomes (SAMN34041181-SAMN34041202) and the metatranscriptomics reads (SRR24029527-SRR24029544) are available in NCBI under BioProject with PRJNA950467 accession. 

Clone repo:
```
$ git clone https://github.com/Chrats-Melkonian/mi_cheese.git
```
See repo structure:
```
$ tree mi_cheese/ -L 3
mi_cheese/
├── code
│   ├── 01_Figure02.R
│   ├── 02_Figure02.R
│   ├── 03_Figure3panelB_E.R
│   ├── 04_community_GEM_sims.sh
│   ├── 04_figure4panelB.R
│   ├── 04_gf_rxns_compile.sh
│   ├── 04_individual_GEM_sims.py
│   ├── 04_model_summary_compile.sh
│   ├── 04_model_summary_generate.py
│   ├── 05_Figure4panelC.R
│   ├── 06_Figure4panelD_E.R
│   ├── 07_Figure5panelA_C.R
│   ├── 08_Figure05PanelE_supp_plots.R
│   ├── 09_figure5panelF.R
│   └── 10_PCAs_supp.R
├── data
│   ├── Models
│   │   ├── curated
│   │   ├── data
│   │   ├── Escher_maps
│   │   └── gapfilled
│   ├── ORF_protein
│   │   ├── LC.faa
│   │   ├── LLm1.faa
│   │   ├── LLm2.faa
│   │   └── ST.faa
│   ├── T0
│   │   ├── Biochemical
│   │   ├── metatranscriptomics
│   │   └── table_analysis
│   └── T1
│       └── Biochemical
├── LICENSE
└── README.md

```

[![DOI](https://zenodo.org/badge/586860014.svg)](https://zenodo.org/badge/latestdoi/586860014)
