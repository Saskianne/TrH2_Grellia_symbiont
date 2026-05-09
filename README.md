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
      - SRR24886407
      - SRR24886411
      - SRR24886417 
      - SRR24886420
      - SRR24886421
      - SRR24886422

### Data Analysis
1. Symbiont gene mapping and counting using CellRanger count pipeline
  - Output: Feature-Barcode-Matrix, Web-Summary  
2. Seurat object handling and data summary
  - 
3. Overlapping with metacell data 

