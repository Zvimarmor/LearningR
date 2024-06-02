
data <- read.table("GSE208783_series_matrix.txt", header = TRUE, fill = TRUE, sep = "\t", quote = "", comment.char = "", skip = 40)

#process the data 
data <- t(data)
identical_columns <- c()

for (i in 1:ncol(data)) {
  current_column <- data[-1, i]
  
  # clean the column
  current_column_processed <- tolower(trimws(current_column))
  actual_data <- current_column_processed[!grepl("^!", current_column_processed) & current_column_processed != ""]
  
  # Check if all values in the column are identical
  if (length(unique(actual_data)) == 1) {
    identical_columns <- c(identical_columns, i)
  }
}

data <- data[, -identical_columns]

print(data)

