# TrH2_Grellia_symbiont
### Code and data from the TrH2 Grellia symbiont analysis  

## Project: "Bioinformatic analysis of single-cell transcriptomics data of Trichoplax sp. H2 on its endosymbiont Grellia incantans"
#### The aim of this project is to find out with which cell type of Trychoplax sp. H2 (TrH2) is associated with the symbionts using single-cell transcriptomics data from Nalje et al. (2023). 

### The provided data were following:
- Symbiont reference genome .fasta, .gtf and .gff files
  - Grellia_incantans_TRH2_HVG.fasta
  - 6666666.192815.gbk.gff
- Single-cell Transcriptomic sequencing data of Trichoplax H2.
  - The scRNA-sequencing data is publicly available in the NCBI data base and can be found by accession numbers. Due to the big files, the data needs to be downloaded with sratoolkit (v3.1.1).
  - scRNA raw data with these (NCBI) accession numbers were analyzed: 
    - SRR24886407: Plac02_H2_H23_10XscRNAseq_10kc
    - SRR24886411:H2H23_4_ACME_10x_10kc
    - SRR24886417: H2drugs_10XscRNAseq_10kc
    - SRR24886420: H2_3_ACME_10x_10kc
    - SRR24886421: H2_2_ACME_10x_10kc
    - SRR24886422: H2_1_ACME_10x_10kc

### Data Analysis
1. Symbiont gene mapping and counting using CellRanger count pipeline
  - Output: Feature-Barcode-Matrix, Web-Summary  
  - Script:
    - First_steps.ipynb
2. Seurat object handling and data summary
  - SEURAT is a tool for single cell genomics analysis, particularly for the clustering analysis.  
  - Scripts: 
    - Features_to_gene.R
    - Seurat_loop_corrected.R
3. Overlapping with metacell data 
  - Metacell clustering (by Najle et al., 2023)
  - Projection of the symbiont signals on the 2D visualization of the single-cell atlas
  - Scripts: 
    - metacell.ipynb
    - xboc_functions.R (from https://github.com/sebepedroslab/Xenoturbella_sc_atlas/)


### Data Files
- TrH2_data_summary_forGit.csv    
  - Summary of the data from metadata from the reference, NCBI, Cellranger count output and seurat analysis
-  absol_gene_count_only.csv
  - Absolute gene counts
- Bact_signal_rel_ab_table_w_average.csv
  - Relative abundance of genes (relative to the total gene count of the library)

### Figures
- See the protocol

