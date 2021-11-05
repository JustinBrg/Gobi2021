library(ggplot2)
library(pcaMethods)
#Read data
zika_table <-read.table(header = FALSE,"Zika.tabular")
count_matrix = zika_table
count_matrix = count_matrix[-1,]
#count_matrix[1] <- NULL
#count_matrix[1:10,1]
#show(count_matrix)

#convert as numeric
count_matrix$V2 <- as.numeric(count_matrix$V2)
count_matrix$V3 <- as.numeric(count_matrix$V3)
count_matrix$V4 <- as.numeric(count_matrix$V4)
count_matrix$V5 <- as.numeric(count_matrix$V5)

head(count_matrix)

#TPM normalization start

#read gene length data
modified_count_matrix <-data.frame(count_matrix)
feature_length <- read.table("featurelength.tabular", sep = "\t", header = T, stringsAsFactors = F)

#divide each read count from the matrix by the given gene length

ordered_feature_length <- feature_length[order(feature_length$Geneid),] #reorder
#head(ordered_feature_length)
ordered_feature_length$Length <- ordered_feature_length$Length/1000
#head(ordered_feature_length)

#modified_count_matrix[1] <- NULL
#head(modified_count_matrix)

number_kb <- ordered_feature_length$Length


modified_count_matrix$V2 <- modified_count_matrix$V2/number_kb
modified_count_matrix$V3 <- modified_count_matrix$V3/number_kb
modified_count_matrix$V4 <- modified_count_matrix$V4/number_kb
modified_count_matrix$V5 <- modified_count_matrix$V5/number_kb
head(modified_count_matrix)

#divided with sum of RPK and 1,000,000 to get TPM
modified_count_matrix$V2 <- modified_count_matrix$V2/(sum(modified_count_matrix$V2)/1E6)
modified_count_matrix$V3 <- modified_count_matrix$V3/(sum(modified_count_matrix$V3)/1E6)
modified_count_matrix$V4 <- modified_count_matrix$V4/(sum(modified_count_matrix$V4)/1E6)
modified_count_matrix$V5 <- modified_count_matrix$V5/(sum(modified_count_matrix$V5)/1E6)


new_new_count_matrix <- modified_count_matrix
head(new_new_count_matrix)



rownames(new_new_count_matrix) <- new_new_count_matrix[,1]
new_new_count_matrix <- new_new_count_matrix[,-1]



#md <- prep(new_new_count_matrix, scale="none", center=TRUE)
pcaa <- pca(t(new_new_count_matrix), scale = "uv",nPcs = 3, method = "svd")

?pca
head(pcaa@network)

head(pcaa@scores)

ggplot(pcaa@scores, aes(PC1, PC2))+ 
  geom_point(size = 3)+
  xlab(paste("PC1 - ", pcaa@R2[1]*100, "%", sep = ""))+
  ylab(paste("PC2 - ",pcaa@R2[2]*100, "%", sep = ""))

#Code from Cheng-Wei Liao
