# MILKOME Cohort Sutdy 
**Date:** 2026.04.27 <br>
**Author:** Yunjeong So

MILKOME is a mother-infant cohort study conducted in Copenhagen, a collaboration between the **Technical University of Denmark (DTU)** and **Rigshospitalet**. The study was led by DTU, with participants recruited by Rigshospitalet.

We sequenced the fecal samples of infants at three different time points: pre-weaning, early-weaning, and late-weaning. For the early-weaning time point, we also sequenced paired fecal samples from the mothers. Additionally, the fecal samples were enriched with diverse glycan sources, and these enriched consortia were sequenced as well.

### Analysis Pipeline

The sequence data analysis pipeline is outlined below. The scripts used can be found under the [`milkome/scripts`](milkome/scripts) directory:

1. Sequence quality check
2. Sequence filtering and host read removal
3. Taxonomic profiling using a *Bifidobacterium longum* subsp. *infantis* and subsp. *longum* resolution database
4. Contig assembly
5. Functional annotation — Carbohydrate Active Enzymes (CAZymes) targeting human milk oligosaccharides (HMOs), dietary fibres (arabinoxylan, pectin, and xyloglucan), and food additives (mannan and gum Arabic)
6. Retrieval of metagenome-assembled genomes (MAGs)
7. Annotation and charactersation of MAGs
8. Downstream analyses

### Repository Structure

```
milkome/
├── scripts/     # Analysis scripts
├── data/        # Metadata and processed data
└── README.md    # Project documentation
```
