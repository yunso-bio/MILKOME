# MILKOME Cohort Sutdy 

MILKOME is a mother-infant cohort study conducted in Copenhagen, a collaboration between the Technical University of Denmark (DTU) and Rigshospitalet. The study was led by DTU, with participants recruited by Rigshospitalet.

We sequenced the fecal samples of infants at three different time points: pre-weaning, early-weaning, and late-weaning. For the early-weaning time point, we also sequenced paired fecal samples from the mothers. Additionally, the fecal samples were enriched with diverse glycan sources, and these enriched consortia were sequenced as well.

The sequence data analysis pipeline is outlined below, and the used scripts can be found under milkome/scripts directory:

1. Sequence quality check
2. Sequence filtering and host read removal
3. Contig assembly
4. Retrieval of metagenome-assembled genomes (MAGs)
5. Mapping MAGs to filtered sequences
6. Annotation and charactersation of MAGs
7. Downstream analyses
