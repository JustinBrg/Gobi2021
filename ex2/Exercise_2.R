library(ggplot2)
library("DESeq2")
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


modified_count_matrix[1] <- NULL



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

#add head
colnames(modified_count_matrix) <-c("Mock", "Mock", "Zika", "Zika")

head(modified_count_matrix)

#clean data remove row with all 0, or it will get error for PCA
remove_all_ziro <-modified_count_matrix[!(modified_count_matrix$Mock== 0 & modified_count_matrix$Zika==0),]

pca <- prcomp(t(remove_all_ziro), scale. = TRUE)
#plot(pca$x[,1],pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)

barplot(pca.var.per, main = "Scree Plot", xlab = "Principal Component", ylab = "Percent Vatiation")

pca.data<- data.frame(Sample = rownames(pca$x), X = pca$x[,1], Y = pca$x[,2])
 head(pca.data)

 #get PCA plot
ggplot(data = pca.data , aes(x = X, y =Y, label = Sample))+
  geom_text(hjust = 0, nudge_x =3)+
  geom_point()+
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = ""))+
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = ""))+
  theme_bw()+
  ggtitle("My PCA Graph")

#try pca Methods in other files


#slplot(pcaa,sl=NULL)



#task 5 DESeq2

new_count_matrix <- count_matrix
head(new_count_matrix)

#new_count_matrix[1] <- NULL 

length(new_count_matrix)

colnames(new_count_matrix) <- c("","SRR3191542.fastq","SRR3191543.fastq","SRR3191544.fastq","SRR3191545.fastq")
rownames(new_count_matrix) <- new_count_matrix[,1]
new_count_matrix <- new_count_matrix[,-1]

#clean low expressed data
row_sub <-apply(new_count_matrix,1,function(row) all(row > 10))
new_count_matrix <- new_count_matrix[row_sub,]

head(new_count_matrix)
head(count_matrix)

expDesign <- data.frame(samples = c("SRR3191542.fastq","SRR3191543.fastq","SRR3191544.fastq","SRR3191545.fastq"), treatment = c("Mock","Mock","Zika","Zika"))


print(expDesign)


dds <-DESeqDataSetFromMatrix( new_count_matrix,
                                 colData = expDesign,
                                 design = ~ treatment)

dds <- DESeq(dds)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="treatment_Zika_vs_Mock")

# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="treatment_Zika_vs_Mock", type="apeglm")



plotMA(res, ylim=c(-2,2))

#Code from Cheng-Wei Liao




