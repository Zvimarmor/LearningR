library(edgeR)
library(ggplot2)
library(ggfortify)
library(biomaRt)
library(matrixStats)

# Read the data
data <- read.table("GSE208783_series_matrix.txt", header = TRUE, fill = TRUE, sep = "\t", skip = 40)

# Process the data 
data <- t(data)
identical_columns <- c()

for (i in seq_len(ncol(data))) {
  current_column <- data[-1, i]
  
  # Clean the column
  current_column_processed <- tolower(trimws(current_column))
  actual_data <- current_column_processed[!grepl("^!", current_column_processed) & current_column_processed != ""]
  
  # Check if all values in the column are identical
  if (length(unique(actual_data)) == 1) {
    identical_columns <- c(identical_columns, i)
  }
}

data <- data[, -identical_columns]
write.csv(data, "GSE208783_series_matrix_processed.csv", row.names = FALSE, quote = FALSE)

# Read the processed series matrix
col_data <- read.csv("GSE208783_series_matrix_processed.csv")

# Read the feature counts data
counts_data <- read.table("GSE208783_gba_feature_counts_rnaseq_april2021.txt", header = TRUE, sep = "\t", quote = "", comment.char = "")

# Remove gene symbol column and transpose counts data
gene_symbols <- counts_data$`Gene symbol`
counts_data <- counts_data[,-1]
rownames(counts_data) <- gene_symbols

# Transpose counts_data
counts_data <- t(counts_data)

# Extract relevant columns from col_data
sample_info <- col_data[2:8, ]  # Adjust rows to match your actual data structure
colnames(sample_info) <- col_data[1, ]  # Use the first row as column names
sample_info <- sample_info[-1, ]  # Remove the first row

# Extract sample IDs from processed data and set them as rownames
sample_info <- as.data.frame(t(sample_info))
sample_ids <- rownames(sample_info)

# Match sample identifiers with counts_data columns
colnames(counts_data) <- sample_ids

# Merge counts_data with sample_info
merged_data <- merge(as.data.frame(counts_data), sample_info, by = "row.names")

# Prepare data for PCA
pca_data <- merged_data[, 2:(ncol(counts_data) + 1)]
pca_data <- as.data.frame(pca_data)
pca_data <- sapply(pca_data, as.numeric)  # Convert to numeric

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Label the PCA plot with 'disease' column
disease_col <- merged_data$disease
autoplot(pca_result, data = merged_data, colour = 'disease_col', label = TRUE, size = 0.3)
autoplot(pca_result, x = 3, y = 4, data = merged_data, colour = 'disease_col', label = TRUE, size = 0.3)

# Prepare data for edgeR analysis
cts <- as.data.frame(counts_data)
col_data <- col_data[-1,]  # Remove the first row with headers
col_data <- subset(col_data, col_data$BUID %in% colnames(cts))

# Order counts and col_data based on sample IDs
cts <- cts[, as.character(col_data$BUID)]
cts <- cts[, order(col_data$BUID)]
col_data <- col_data[order(col_data$BUID), ]

# Check rownames
stopifnot(rownames(col_data) == colnames(cts))

# Create DGEList object
y <- DGEList(counts = cts, group = col_data$condition)

# Filter genes by expression
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

# Create normalized counts matrix
cts_norm <- as.data.frame(cpm(y, log = FALSE))

# Model matrix for differential expression analysis
dsgn <- model.matrix(~ RIN + batch + sex + condition, data = col_data)
y <- estimateDisp(y, dsgn, robust = TRUE)

# Fit the model and perform likelihood ratio test
fit <- glmQLFit(y, dsgn, robust = TRUE)
lrt <- glmLRT(fit, coef = 5)

# Extract top differentially expressed genes
top_genes <- as.data.frame(topTags(lrt, adjust.method = 'fdr', n = nrow(cts_norm)))
top_genes$transcript <- rownames(top_genes)

# Print the top differentially expressed genes
print(head(top_genes))
