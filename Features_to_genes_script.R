## SEURAT object feature ID to gene name in the Annoatation data

# Load necessary library
library(dplyr)
library(stringr)

# Read the file line by line
gff_file <- "D:/scTriH2/symbiont_ref_genome/6666666.192815.gbk.gff"
gff_lines <- readLines(gff_file)

# Find where the FASTA section starts (if it exists)
fasta_start <- grep("^##FASTA", gff_lines)

# If FASTA exists, exclude lines after it; otherwise, read the whole file
if (length(fasta_start) > 0) {
  gff_lines <- gff_lines[1:(fasta_start - 1)]
}

# Write the GFF part to a temporary file for reading as a table
temp_gff_file <- tempfile(fileext = ".gff")
writeLines(gff_lines, temp_gff_file)

# Read the GFF content
gff_data <- read.table(temp_gff_file, sep = "\t", header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE)

# Assign GFF column names
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Check the first few rows to ensure it's read correctly
head(gff_data)

# Filter for CDS features
cds_data <- gff_data %>% filter(feature == "CDS")

# Extract IDs and Products
cds_data <- cds_data %>%
  mutate(
    ID = str_extract(attribute, "ID=([^;]+)"),          # Extract ID
    product = str_extract(attribute, "product=([^;]+)") # Extract product
  ) %>%
  mutate(
    ID = sub("ID=", "", ID),                            # Clean ID
    product = sub("product=", "", product)              # Clean product name
  )

# Check the mapping
head(cds_data[, c("ID", "product")])



# Mapping IDs to Seurat object features
id_to_gene <- setNames(cds_data$product, cds_data$ID)

# Adjust the feature IDs in the mapping to match Seurat's format (replace underscores with dashes)
id_to_gene_adjusted <- id_to_gene
names(id_to_gene_adjusted) <- gsub("_", "-", names(id_to_gene))

# Extract current feature names from the merged Seurat object
current_features <- rownames(SRR24886410_RDS)

# Initialize updated features vector
updated_features <- current_features

# Loop to replace features with real gene names if they exist in the adjusted mapping
for (i in seq_along(current_features)) {
  feature_id <- current_features[i]
  
  # Replace with gene name if there's a match
  if (feature_id %in% names(id_to_gene_adjusted)) {
    updated_features[i] <- id_to_gene_adjusted[[feature_id]]
  } else {
    updated_features[i] <- feature_id  # Keep original if no match
  }
}

# Update the rownames in the merged Seurat object
rownames(SRR24886410_RDS) <- updated_features

# Verify changes
head(rownames(SRR24886410_RDS))



