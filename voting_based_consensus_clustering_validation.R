# Freie Universitaet Berlin
# Fachbereich Mathematik und Informatik
# Institut fuer Informatik, Institut fuer Bioinformatik

# Calvinna Yosephine Caswara

#------------------Packages---------------
library(dendextend)
library(ape)
library(dendextend)
library(colorspace)
library(gplots)
library(corrplot)
unloadNamespace("dendextend"); attachNamespace("dendextend")

#------------------Initialization---------------
pathrunscript = "C:/Users/hp/ownCloud/Documents/Software/"
mainDir <- "C:/Users/hp/ownCloud/Documents/Software/"

# create output folder
subDir <- "output/"
dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))
pathrda = "C:/Users/hp/ownCloud/Documents/Software/output/"

# create output subfolders
subDir_k <- "determinedk/"
subDir_cc <- "consensusclustering/"
subDir_gep <- "geneexprprofiles/"

dataset="glioblastoma"
linkage="complete"           # for the dendrogram
mng=10

#------------------Packages and functions--------------- 
source(paste(pathrunscript,"newwave_takeout-notAdjRand-versionRMP-function.R",sep=""))
data(pwdata.pathways.KEGG.hsa)

source(paste(pathrunscript,"voting_based_consensus_clustering_functions.R",sep=""))

#------------------read expression data and annotation of expression data--------------- 
#Loading glioblastoma Data
exprdata1 <- read.table("~/Bachelorarbeit/Software/dataset/gb-exprdata-eumrc-norm_DESeq_1step-pseudocount-log2-asEntrez-unique-max-sorted.tsv",header=T, row.names=1)
glioblastoma <- read.table("~/Bachelorarbeit/Software/dataset/samples2class-all-sorted.tsv",header=T, row.names=1, fill=TRUE)
#Loading neuroblastoma Data
exprdata2 <- read.table("~/Bachelorarbeit/Software/dataset/neuroblastoma-MAQC-1col/MAavg_modified-matrix-asEntrez-viaEnsembl70.tsv",header=T, row.names=1, fill=TRUE)
neuroblastoma <- read.table("~/Bachelorarbeit/Software/neuroblastoma-MAQC-1col/Cologne.col2and13.tsv",header=T, fill=TRUE)
#Loading pediatric tumor Data
exprdata3 <- read.table("~/Bachelorarbeit/Software/dataset/brain_tumors_pediatric-GSE35493/brain_tumors_pediatric-GSE35493/GSE35493_series_matrix-cleaned-uniqueEntrez-sorted.tsv",header=T, row.names=1, fill=TRUE)
pediatric <- read.table("~/Bachelorarbeit/Software/brain_tumors_pediatric-GSE35493/brain_tumors_pediatric-GSE35493/GSE35493_series_matrix-samples2clinical_annotation-sorted.tsv",header=TRUE, fill=TRUE, sep="\t")
#Loading breast cancer Data
exprdata4 <- read.table("~/Bachelorarbeit/Software/dataset/Breast_cancer-GSE48390-various-Taipeh-81_tumors/Breast_cancer-GSE48390-various-Taipeh-81_tumors/GSE48390_series_matrix-cleaned-uniqueEntrez.tsv",header=T, row.names=1, fill=TRUE)
breastcancer <- read.table("C:/Users/hp/ownCloud/Documents/Software/dataset/Breast_cancer-GSE48390-various-Taipeh-81_tumors/Breast_cancer-GSE48390-various-Taipeh-81_tumors/annotation_nbs.txt",header=TRUE, fill=TRUE, sep="\t")
#Loading neuroblastoma primary tumors Data
exprdata5 <- read.table("~/Bachelorarbeit/Software/dataset/Neuroblastoma-GSE16237-51_mostly_primary_tumors/Neuroblastoma-GSE16237-51_mostly_primary_tumors/GSE16237_series_matrix-cleaned-uniqueEntrez.tsv",header=T, row.names=1, fill=TRUE)
neuroblastomaprimary <- read.table("C:/Users/hp/ownCloud/Documents/Software/dataset/Neuroblastoma-GSE16237-51_mostly_primary_tumors/Neuroblastoma-GSE16237-51_mostly_primary_tumors/annotation_nbs.txt",header=T, fill=TRUE)


#------------------prepare annotation for comparison--------------- 
if(dataset=="glioblastoma"){
  actual_nb_clusterings = 3
  exprdata = exprdata1
}else if(dataset=="neuroblastoma"){
  actual_nb_clusterings = 7
  exprdata = exprdata2
}else if(dataset=="pediatric"){
  actual_nb_clusterings = 5
  exprdata = exprdata3
}else if(dataset=="breastcancer"){
  actual_nb_clusterings = 4
  exprdata = exprdata4
}else if(dataset=="neuroblastomaprimary"){
  actual_nb_clusterings = 5
  exprdata = exprdata5
}

ktest=actual_nb_clusterings

  
if(dataset=="glioblastoma"){
  consensusclustering_original <- as.character(1:3)[ match(glioblastoma[,1], c("1_wt","2_mut","3_control"))]
  consensusclustering_original <- c(as.integer(consensusclustering_original))
  names(consensusclustering_original) <- rownames(glioblastoma)
  tumor_type = c("IDH wildtype","IDH1 mutation","normal brain")
  
}else if(dataset=="neuroblastoma"){
  consensusclustering_original <- gsub('2A', '5', neuroblastoma[,2])
  consensusclustering_original <- gsub('2B', '6', consensusclustering_original)
  consensusclustering_original <- gsub('4S', '7', consensusclustering_original)
  consensusclustering_original <- gsub('localised_multilocular', '1', consensusclustering_original)
  consensusclustering_original <- as.integer(consensusclustering_original)
  names(consensusclustering_original) <- neuroblastoma[,1]
  digits = paste("X", names(consensusclustering_original[grep(pattern="^[[:digit:]]", names(consensusclustering_original))]), sep="")
  nodigits = names(consensusclustering_original[grep(pattern="[^[[:digit:]]]*", names(consensusclustering_original))])
  names(consensusclustering_original)=c(digits,nodigits)
  tumor_type = c("stage 1","stage 2","stage 3","stage 4","stage 2A","stage 2B","stage 4S")
  
}else if(dataset=="pediatric"){
  consensusclustering_original <- gsub('primitive neuroectodermal tumor', '1', pediatric[,2])
  consensusclustering_original <- gsub('normal brain', '2', consensusclustering_original)
  consensusclustering_original <- gsub('normal frontal', '2', consensusclustering_original)
  consensusclustering_original <- gsub('normal cerebellum', '2', consensusclustering_original)
  consensusclustering_original <- gsub('medulloblastoma, large-cell', '3', consensusclustering_original)
  consensusclustering_original <- gsub('medulloblastoma', '3', consensusclustering_original)
  consensusclustering_original <- gsub('glioblastoma', '4', consensusclustering_original)
  consensusclustering_original <- gsub('atypical teratoid / rhabdoid tumor', '5', consensusclustering_original)
  consensusclustering_original <- as.integer(consensusclustering_original)
  names(consensusclustering_original) <- pediatric[,1]
  tumor_type = c("primitive neuroectodermal tumor","normal brain","medulloblastoma","glioblastoma","ATRT" )
  
}else if(dataset=="breastcancer"){
  consensusclustering_original <- breastcancer[,2]
  names(consensusclustering_original) <- breastcancer[,1]
  tumor_type = c("ER-positive","normal breast-like","HER2-positive","HER2- and ER-pos.")
  
}else if(dataset=="neuroblastomaprimary"){
  consensusclustering_original <- gsub('4S', '5', neuroblastomaprimary[,2])
  consensusclustering_original <- as.integer(consensusclustering_original)
  names(consensusclustering_original) <- neuroblastomaprimary[,1]
  
  tumor_type = c("stage 1","stage 2","stage 3","stage 4","stage 4S")
}

###########################################################
# 
# VALIDATION
# 
###########################################################


load(paste(pathrda, subDir_cc ,dataset, "_consensusclustering_euclid_k",ktest,"_minpw",mng,".rda", sep =""))
load(paste(pathrda, subDir_cc ,dataset, "_consensusclustering_spearman_k",ktest,"_minpw",mng,".rda", sep =""))
load(paste(pathrda, subDir_cc ,dataset, "_consensusclustering_pearson_k",ktest,"_minpw",mng,".rda", sep =""))

valpath=paste(pathrda, subDir_cc, sep="")

#------------------RI/ARI/NMI Table from k=2...15--------------- 

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
distance ="euclid"
corrplot(as.matrix(valmatagg(valpath, dataset, distance,consensusclustering_original,mng)), method="color", col=col(200), 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=90, number.digits =3 #Text label color and rotation
)
distance ="spearman"
corrplot(as.matrix(valmatagg(valpath, dataset, distance,consensusclustering_original,mng)), method="color", col=col(200), 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=90, number.digits =3 #Text label color and rotation
)
distance ="pearson"
corrplot(as.matrix(valmatagg(valpath, dataset, distance,consensusclustering_original,mng)), method="color", col=col(200), 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=90, number.digits =3 #Text label color and rotation
)

#------------------Levelplot: Consensus Clusterin on levels/labels--------------- 

x=1:length(consensusclustering_original)
if(dataset=="glioblastoma"){
  Colour <- cut(x, breaks = c(0, 8,16,24),labels = rainbow_hcl(length(tumor_type)))
}else if(dataset=="neuroblastoma"){
  Colour <- cut(x, breaks = c(0, 119,122,193,341,382,416,478),labels = rainbow_hcl(length(tumor_type)))
}else if(dataset=="pediatric"){
  Colour <- cut(x, breaks = c(0, 9,18,39,51,71),labels = rainbow_hcl(length(tumor_type)))
}else if(dataset=="breastcancer"){
  Colour <- cut(x, breaks = c(0, 37,47,63,81),labels = rainbow_hcl(length(tumor_type)))
}else if(dataset=="neuroblastomaprimary"){
  Colour <- cut(x, breaks = c(0, 21,27,32,45,51),labels = rainbow_hcl(length(tumor_type)))
}

consensusclustering_original[order(consensusclustering_original)]


par(mar = c(5,6,2,14))
plot(consensusclustering_euclid$clustering[order(consensusclustering_original)], col=as.character(Colour), xlab="Sample", ylab="Cluster Label",pch=19) 
par(xpd=TRUE)
op <- par(cex=.64)
legend("right",inset=c(-0.2,-0.5), horiz=FALSE, legend = tumor_type, col = rainbow_hcl(length(tumor_type)),pch=19,cex = 1.3)
for(type in 1:(length(tumor_type)-1)){
  mx = max(which(consensusclustering_original[order(consensusclustering_original)]==type))
  clip(1,length(consensusclustering_original[order(consensusclustering_original)]),1,length(tumor_type))
  abline(v=mx+0.5)
}
par(op)

par(mar = c(5,6,2,14))
plot(consensusclustering_spearman$clustering[order(consensusclustering_original)], col = as.character(Colour), xlab="Sample", ylab="Cluster Label",pch=19)  
par(xpd=TRUE)
op <- par(cex=.64)
legend("right",inset=c(-0.2,-0.5), horiz=FALSE, legend = tumor_type, col = rainbow_hcl(length(tumor_type)),pch=19,cex = 1.3)
for(type in 1:(length(tumor_type)-1)){
  mx = max(which(consensusclustering_original[order(consensusclustering_original)]==type))
  clip(1,length(consensusclustering_original[order(consensusclustering_original)]),1,length(tumor_type))
  abline(v=mx+0.5)
}
par(op)

par(mar = c(5,6,2,14))
plot(consensusclustering_pearson$clustering[order(consensusclustering_original)], col = as.character(Colour), xlab="Sample", ylab="Cluster Label",pch=19)  
par(xpd=TRUE)
op <- par(cex=.64)
legend("right",inset=c(-0.2,-0.5), horiz=FALSE, legend = tumor_type, col = rainbow_hcl(length(tumor_type)),pch=19,cex = 1.3)
for(type in 1:(length(tumor_type)-1)){
  mx = max(which(consensusclustering_original[order(consensusclustering_original)]==type))
  clip(1,length(consensusclustering_original[order(consensusclustering_original)]),1,length(tumor_type))
  abline(v=mx+0.5)
}
par(op)

#------------------Dendrograms/Hierarchical Clustering--------------- 

valpath2=paste(pathrda, subDir_gep, sep="")
load(paste(valpath2, dataset, "_pw_list_copy_filt_minpw",mng,".rda", sep =""))
load(paste(valpath2, dataset, "_clustering_height_minpw",mng,".rda", sep =""))

#------------------Euclidean Distance---------------
dista="euclid"
load(paste(valpath, dataset,"_votingmatrixe_k",actual_nb_clusterings,"_minpw",mng,".rda", sep =""))
dissimilaritymatrix <- as.dist(1-similaritymatrixe/clustering_height)
similaritymatrix = similaritymatrixe
hc_cons <- hclust(dissimilaritymatrix, method = linkage)
consensusclustering_original_lev <- unique(consensusclustering_original)

dend <- as.dendrogram(hc_cons)
dend <- rotate(dend, 1:length(consensusclustering_original))
tumtypeorder=1:length(tumor_type)

labels_colors(dend) <-  rainbow_hcl(actual_nb_clusterings)[consensusclustering_original[order.dendrogram(dend)]]
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.5)
par(mar = c(10,3,3,7))
plot(dend, main = "", horiz =  TRUE,  nodePar = list(cex = .007))
par(xpd=TRUE)
op <- par(cex=.64)
legend("bottom", inset=c(0,-0.25), legend = tumor_type, fill = rainbow_hcl(actual_nb_clusterings))
rect.dendrogram(dend, k=actual_nb_clusterings, horiz = TRUE, border = 8, lty = 5, lwd = 2)
par(op)
#------------------Spearman Distance---------------
dista="spearman"
load(paste(valpath, dataset,"_votingmatrixs_k",actual_nb_clusterings,"_minpw",mng,".rda", sep =""))
dissimilaritymatrix <- as.dist(1-similaritymatrixs/clustering_height)
similaritymatrix = similaritymatrixs
hc_cons <- hclust(dissimilaritymatrix, method = linkage)
consensusclustering_original_lev <- unique(consensusclustering_original)

dend <- as.dendrogram(hc_cons)
dend <- rotate(dend, 1:length(consensusclustering_original))
tumtypeorder=1:length(tumor_type)

labels_colors(dend) <-  rainbow_hcl(actual_nb_clusterings)[consensusclustering_original[order.dendrogram(dend)]]
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.5)
par(mar = c(10,3,3,7))
plot(dend, main = "", horiz =  TRUE,  nodePar = list(cex = .007))
par(xpd=TRUE)
op <- par(cex=.64)
legend("bottom", inset=c(0,-0.25), legend = tumor_type, fill = rainbow_hcl(actual_nb_clusterings))
rect.dendrogram(dend, k=actual_nb_clusterings, horiz = TRUE, border = 8, lty = 5, lwd = 2)
par(op)
#------------------Pearson Distance---------------
dista="pearson"
load(paste(valpath, dataset,"_votingmatrixp_k",actual_nb_clusterings,"_minpw",mng,".rda", sep =""))
dissimilaritymatrix<- as.dist(1-similaritymatrixp/clustering_height)
similaritymatrix = similaritymatrixp
hc_cons <- hclust(dissimilaritymatrix, method = linkage)
consensusclustering_original_lev <- unique(consensusclustering_original)

dend <- as.dendrogram(hc_cons)
dend <- rotate(dend, 1:length(consensusclustering_original))
tumtypeorder=1:length(tumor_type)

labels_colors(dend) <-  rainbow_hcl(actual_nb_clusterings)[consensusclustering_original[order.dendrogram(dend)]]
dend <- hang.dendrogram(dend,hang_height=0.1)
dend <- set(dend, "labels_cex", 0.5)
par(mar = c(10,3,3,7))
plot(dend, main = "", horiz =  TRUE,  nodePar = list(cex = .007))
par(xpd=TRUE)
op <- par(cex=.64)
legend("bottom", inset=c(0,-0.25), legend = tumor_type, fill = rainbow_hcl(actual_nb_clusterings))
rect.dendrogram(dend, k=actual_nb_clusterings, horiz = TRUE, border = 8, lty = 5, lwd = 2)
par(op)

#------------------Heatmap of the Distance Matrix--------------- 
corrplot(as.matrix(dissimilaritymatrix)[order(consensusclustering_original),order(consensusclustering_original)], method="color", tl.cex=0.5, tl.col="black")


