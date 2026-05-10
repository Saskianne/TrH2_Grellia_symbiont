library(svMisc)
library(dplyr)
library(Seurat)
library(patchwork)
library(BiocManager)
library(Biostrings)
library(Matrix)
library(stringr)
library(rtracklayer)


parent_directory <- "D:/scTriH2/SEURAT/extracted/1.cellRangerCount"
setwd(parent_directory)

options(Seurat.object.assay.version = 'v4')

# # Assign GFF column names
GFFcolnames()
  # c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff_path <- "D:/scTriH2/symbiont_ref_genome/6666666.192815.gbk.gff"
gff_data <- readGFF(gff_path)
gff_table <- readGFF(gff_path, tags = character(0))

head(gff_data)
head(gff_table)

# Filter for features of interest: CDS, rRNA, tRNA
CDS_data <- gff_table %>% filter(type == "CDS")
rRNA_data <- gff_table %>% filter(type == "rRNA")
tRNA_data <- gff_table %>% filter(type == "tRNA")
exon_data <- gff_table%>% filter(type == "exon")
feature_data_list <- list(CDS_data, rRNA_data, tRNA_data)
id_to_transcript <- c()
hypothetical_protein_counter <- 1

for (feature_data in feature_data_list){
  feature_data <- feature_data %>%
    # Extract ID field and gene name 
    mutate(
      ID = str_extract(attributes, "ID=([^;]+)"),                   
      product = str_extract(attributes, "product=([^;]+)")         
    ) %>%
    # clean the ID and gene name
    mutate(
      ID = sub("ID=", "", ID),                                    
      product = sub("product=", "", product)                      
    ) %>%
    # leave the name features if product is an hypothetical protein
    mutate(
      product = ifelse(product == "hypothetical protein", ID, product)
    ) %>%
    # Add suffix if product name repeats
    mutate(
      product = ave(product, product, FUN = function(x) {
        if (length(x) > 1) paste0(x, "_", seq_along(x)) else x
      })
    )

  
  
  head(feature_data[, c("ID", "product")])
  id_to_transcript <- c(id_to_transcript, setNames(feature_data$product, feature_data$ID))
} 

id_to_transcript <- unlist(id_to_transcript)
head(id_to_transcript)

# List of SRRs to work with
SRR_list_full <- list.dirs(parent_directory, full.names = FALSE, recursive = FALSE)
SRR_list_full <- SRR_list_full[grepl("SRR*", SRR_list_full)]
# exclude the SRR nr. without any bacterial gene reads
SRR_list <- SRR_list_full[!SRR_list_full %in% c("SRR24886387", "SRR24886410")] #, "SRR24886416"
cat(SRR_list)
sum_list <- list()
seurat_list <- list()
so_list <- list()
list_for_tab <- c(SRR_list, sum_list)

# Output directory
# dir.create("Seurat_outputs", showWarnings = FALSE)
base_dir <- file.path(parent_directory, "Seurat_outputs")
dir.create(base_dir, showWarnings = FALSE)

# Loop for filtering and seurat object starts
cat("Loop starts \n")
for (SRR in SRR_list) {
  cat("Working with", SRR, "\n")
  SRR_output <- paste0(SRR, "_output")
  output_dir <- file.path(base_dir, SRR_output)
  
  # Read 10X data
  raw_dir <- file.path(parent_directory, SRR, "outs/raw_feature_bc_matrix/")
  raw_matrix <- Read10X(data.dir = raw_dir)
  cat("Created the matrix: raw_matrix_", SRR, "\n", file = outfile_conn)
  
  # Data filtering
  # Filtering by min. mapped read count per cell
  raw_matrix_filtered <- raw_matrix[rowSums(raw_matrix > 0) >= 1, ]
  # Extract current feature names from the merged Seurat object
  current_features <- rownames(raw_matrix_filtered)
  # Initialize updated features vector
  updated_features <- current_features
  
  # Replace row names in Seurat object with gene names
  updated_features <- ifelse(current_features %in% names(id_to_transcript), 
                                     id_to_transcript[current_features], 
                                     current_features)
  
  # Update the rownames in the merged Seurat object
  rownames(raw_matrix_filtered) <- updated_features
  
  head(rownames(raw_matrix_filtered))
  
  cat("Chekpoint1")
  # change the names of the barcode names according to the sample 
  if (SRR == "SRR24886407") {
    colnames(raw_matrix_filtered) <- paste0("Plac02_H2_H23_10XscRNA_", colnames(raw_matrix_filtered))
  } else if (SRR == "SRR24886411") {
    colnames(raw_matrix_filtered) <- paste0("H2H23_4_ACME_10x_10kc_", colnames(raw_matrix_filtered))
  } else if (SRR == "SRR24886416" | SRR == "SRR24886417") {
    colnames(raw_matrix_filtered) <- paste0("H2drugs_10XscRNAseq_10kc_", colnames(raw_matrix_filtered))
  } else if (SRR == "SRR24886420") {
    colnames(raw_matrix_filtered) <- paste0("H2_3_ACME_10x_10kc_", colnames(raw_matrix_filtered))
  } else if (SRR == "SRR24886421") {
    colnames(raw_matrix_filtered) <- paste0("H2_2_ACME_10x_10kc_", colnames(raw_matrix_filtered))
  } else if (SRR == "SRR24886422") {
    colnames(raw_matrix_filtered) <- paste0("H2_1_ACME_10x_10kc_", colnames(raw_matrix_filtered))
  }
  
  gene_rank <- data.frame(
    Gene = rownames(raw_matrix_filtered),
    Counts = rowSums(raw_matrix_filtered)
  )
  
  # rank genes by their abundance
  gene_names <- rownames(raw_matrix)
  # gene_rank_all <- data.frame(Gene = gene_names, SRR = RowSum)
  gene_ranking <- gene_rank[order(gene_rank$Counts, decreasing = TRUE), ]
  high_genes <- head(gene_ranking, 50)
  write.table(high_genes, file = file.path(output_dir, "high_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  low_genes <- gene_ranking[gene_ranking$Counts <="3",]
  write.table(low_genes, file = file.path(output_dir, "low_genes.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)  
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(counts = raw_matrix_filtered, project = SRR, 
                                      min.features = 1, min.cells = 1)
  # make a seurat list for creating a gene count table
  assign(paste0("object_", SRR), seurat_object)
  
  i <- match(SRR, SRR_list)
  so_list <- append(so_list, paste0("object_", SRR), after = i-1 )
  seurat_list <- append(seurat_list, assign(paste0("object_", SRR), seurat_object), after=i-1 )
  
  # Add AvReadsPerGenes(average coverage of a gene in a cell) as layer 
  # seurat_object$GenesPerCell <- (seurat_object$nFeature_RNA / seurat_object$nCount_RNA)
  seurat_object$AvReadsPerGenes <- (seurat_object$nCount_RNA / seurat_object$nFeature_RNA)
  
  # sum all bacterial counts in the sample 
  sum_counts <- sum(rowSums(seurat_object[["RNA"]]$counts))
  i <- match(SRR, SRR_list)
  
  sum_list <- append(sum_list, sum_counts, after= i-1)
  cat("The total read count is ", sum_counts , "\n", file = outfile_conn)
  
  # # gene_rank of seurat_object
  # so_gene_rank <- data.frame(
  #   Gene = rownames(seurat_object[["RNA"]]$counts),
  #   RowSum = rowSums(seurat_object[["RNA"]]$counts)
  # )
  # so_gene_ranking <- so_gene_rank[order(so_gene_rank$RowSum, decreasing = TRUE), ]
  # so_high_genes <- head(so_gene_ranking, 50)
  # write.table(so_high_genes, file = file.path(output_dir, "SO_high_genes.tsv"),
  #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  # so_low_genes <- so_gene_rank[so_gene_rank$RowSum <="3",]
  # write.table(so_low_genes, file = file.path(output_dir, "SO_low_genes.tsv"),
  #             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE) 
  
  # Save filtered matrix in 10X format (in 3 tsv files)
  writeMM(GetAssayData(seurat_object, layer = "counts"), file = file.path(output_dir, "matrix.mtx"))
  
  features_data <- data.frame(
    Feature_ID = rownames(seurat_object),
    Feature_Name = rownames(seurat_object),  
    Feature_Type = "Gene Expression"
  )
  write.table(features_data, file = file.path(output_dir, "features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(seurat_object), file = file.path(output_dir, "barcodes.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Gzip the files to match Read10X expectations
  files_to_gzip <- c("matrix.mtx", "features.tsv", "barcodes.tsv")
  for (file in files_to_gzip) {
    gz_file <- gzfile(file.path(output_dir, paste0(file, ".gz")), "w")
    writeLines(readLines(file.path(output_dir, file)), gz_file)
    close(gz_file)
    file.remove(file.path(output_dir, file))  # Remove the original uncompressed file
  }
  
  cat("The filtered matrix is saved in 10X gzipped format.\n")
  
  # Save as RDS file
  RDS_name <- paste0(SRR, "_filtered.rds")
  SaveSeuratRds(seurat_object, file = file.path(output_dir,RDS_name))
}


# Create a table of gene counts
names(seurat_list) <- so_list
names(sum_list) <- so_list
gene_table <- lapply(seurat_list, function(obj) {
  counts <- rowSums(obj[["RNA"]]$counts)
  data.frame(Gene = rownames(obj[["RNA"]]$counts), Counts = counts)
})
# Combine into one data frame by gene name
combined_table <- Reduce(function(x, y) full_join(x, y, by = "Gene"), gene_table)
# Replace NAs with 0 (in case of missing genes in some samples)
combined_table[is.na(combined_table)] <- 0
colnames(combined_table) <- c("Gene", names(seurat_list))
# Export to CSV
write.csv(combined_table, "gene_counts_table.csv", row.names = FALSE)
# Calculate sum_counts for each sample
sum_counts <- sapply(seurat_list, function(obj) sum(rowSums(obj[["RNA"]]$counts)))
# Normalize counts so all samples have the same total counts
target_sum <- 30000 
normalized_table <- combined_table %>%
  mutate(across(-Gene, ~ . / sum_counts[which(names(sum_counts) == cur_column())] * target_sum))
# Export normalized counts
write.csv(normalized_table, "normalized_gene_counts.csv", row.names = FALSE)

# Rank genes by total counts across samples
ranked_table <- normalized_table %>%
  mutate(TotalCounts = rowSums(across(-Gene))) %>%
  arrange(desc(TotalCounts))
# Export ranked table
write.csv(ranked_table, "ranked_normalized_gene_counts.csv", row.names = FALSE)

read

# exclude the 4th column with very low count
file_path <- "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/ranked_normalized_gene_counts.csv"
gene_data <- read_csv(file_path)

gene_data_filtered <- gene_data[, c(-4,-9)]
sample_cols <- 2:(ncol(gene_data_filtered) - 1)
gene_data_filtered$Sum <- rowSums(gene_data_filtered[, sample_cols])
gene_data_filtered <- gene_data_filtered %>% arrange(desc(Sum))

write_csv(gene_data_filtered, "/mnt/d/scTriH2/SEURAT/extracted/1.cellRangerCount/ranked_norm_gene_counts_high_only.csv")


cat("Run for every SRR in SRR_list is done \n")

# Create a table with bacterial and host 
# Combine lists into a data frame
df <- data.frame(
  SRR = unlist(SRR_list),
  SUM_counts = unlist(sum_list),
)
write.csv(df, "abundance_table.csv", row.names = FALSE)



cat("Start working with merging \n")


# Assuming SRR_paths contains the correct paths to each SRR output directory.
# Ensure SRR_paths is correctly defined:
SRR_paths <- file.path("D:/scTriH2/SEURAT/extracted/1.cellRangerCount/Seurat_outputs", paste0(SRR_list, "_output"))

# Check if all directories exist before proceeding
if (!all(dir.exists(SRR_paths))) {
  stop("Some SRR paths do not exist. Please verify the paths in SRR_paths.")
}

# Read matrices and create Seurat objects
seurat_list <- lapply(SRR_paths, function(path) {
  mat <- Read10X(data.dir = path)
  CreateSeuratObject(counts = mat, project = basename(path), min.cells = 2, min.features = 2)
})

for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$Run <- SRR_list[i]  # Add SRR ID as metadata
}

# Merge Seurat objects
summary_matrix_1 <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1], 
  collapse = TRUE,
  Project = getOption(x = "scTrH2_SP_1")
)

cat("Merging done \n")

# Save the merged object
saveRDS(summary_matrix, file = "E:/scTriH2/SEURAT/extracted/1.cellRangerCount/Seurat_outputs/merged_seurat_object.rds")
cat("Merged Seurat object saved successfully.\n")



