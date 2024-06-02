library(edgeR) 
library(ggplot2)
library(ggfortify)
library('biomaRt')
library(matrixStats)


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
  if (length(unique(actual_data)) == 1){
    identical_columns <- c(identical_columns, i)
  }
}

data <- data[, -identical_columns]

# Save the processed data and convert it to csv format
write.csv(data, "GSE208783_series_matrix_processed.csv", row.names = FALSE)

# Define CountsDataFrame based on processed data
CountsDataFrame <- data
ColData <- read.table("GSE208783_series_matrix_processed.csv", header = TRUE, sep = ",", quote = "", comment.char = "", skip = 0)

# Check for duplicated Sample_ID in pcaData1
duplicated_ids <- pcaData1$Sample_ID[duplicated(pcaData1$Sample_ID)]
print(duplicated_ids)

tmp1<-CountsDataFrame
pcaData1<-as.data.frame(t(tmp1))
pcaData1$Sample_ID<-rownames(pcaData1)
pcaData1<-merge(pcaData1,ColData,by='Sample_ID')
rownames(pcaData1)<-pcaData1$Sample_ID 
pca_res1<-prcomp(pcaData1[,2:nrow(tmp1)],scale. = T)

c='RIN' # change c to the parameter you want to lable the PCA with
autoplot(pca_res1, x=1,y=2 , data = pcaData1, colour = c,label=T,size=0.3)
autoplot(pca_res1, x=3,y=4 , data = pcaData1, colour = c,label=T,size=0.3)


cts<-CountsDataFrame
cold1<-ColData
cold1<-subset(cold1,cold1$BUID %in% colnames(cts))
cts<-cts[,as.character(cold1$BUID)] 
cts<-cts[,order(cold1$BUID)]
cold1<-cold1[order(cold1$BUID),]
row(cold1)==sum(cold1$BUID==colnames(cts))
y <- DGEList(counts=cts,group=cold1$condition) # change the condition to the name of your condition column
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

cts1<-as.data.frame(cpm(y,log = F)) # creates the normalized counts matrix

dsgn <- model.matrix(~RIN+batch+sex+condition, data = cold1)
y <- estimateDisp(y, dsgn, robust = T)

## if you want, you can run the next code to see if the different coeficients correlate with eachother:
# logFC <- predFC(y,dsgn,prior.count=1,dispersion=0.05) ; cor(logFC) ; plotBCV(y)

head(dsgn)
# you change the coef to the coeficient that intersts you (the number of the column in the dsgn matrix)
fit <- glmQLFit(y, dsgn, robust = T)
lrt1 <- glmLRT(fit,coef = 5) 

sgGens<-as.data.frame(topTags(lrt,adjust.method = 'fdr',n = nrow(cts1)))
sgGens$transcript<-rownames(sgGens)



