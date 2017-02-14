# Freie Universitaet Berlin
# Fachbereich Mathematik und Informatik
# Institut fuer Informatik, Institut fuer Bioinformatik

# Calvinna Yosephine Caswara

#------------------Packages and functions--------------- 
pathrunscript = "C:/Users/hp/ownCloud/Documents/Software/"
source(paste(pathrunscript,"newwave_takeout-notAdjRand-versionRMP-function.R",sep=""))
data(pwdata.pathways.KEGG.hsa)

source(paste(pathrunscript,"voting_based_consensus_clustering_functions.R",sep=""))

library(factoextra)
library(fpc)
library(cluster)
library(NbClust)
library(PathWave)
library(plyr)
library(bioDist)
library(clusterCrit)
library(flexclust)
library(clue)
library(ape)
library(dendextend)
library(colorspace)

#------------------Initialization--------------- 
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
dir.create(file.path(pathrda, subDir_k))
setwd(file.path(pathrda, subDir_k))
dir.create(file.path(pathrda, subDir_cc))
setwd(file.path(pathrda, subDir_cc))
dir.create(file.path(pathrda, subDir_gep))
setwd(file.path(pathrda, subDir_gep))


dataset="glioblastoma"
mng=10
nbmax=15

#------------------Rename Data---------------
#imported datasets from KEGG
pw_id_table <- pwdata.pid2pname.KEGG.hsa
pw_gene_table <- pwdata.reac.genes.KEGG.hsa

dim_pw_id <- dim(pw_id_table)
dim_pw_gene <- dim(pw_gene_table)

list_of_pw <- vector("list", dim_pw_id[1]) 
table_pw_genes <- vector("list", dim_pw_id[1]) 

elist <- c("hsa01100", strsplit(paste("hsa0",2000:5416, sep=""), " "))
indices_clustering_validation <- c("kl", "ch", "hartigan", "cindex", "db",  "duda", "ratkowsky", "ptbiserial", "gamma", "tau", "dunn", "sdindex")



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


###########################################################
# 
# MAIN PART
# 
###########################################################

#------------------Preprocessing Data----------------

print("preprocessing data")

#exclude too big or irrelevant pathways first
pw_id_table_filtered <- excludepw(pw_id_table, elist, print=FALSE)


#reads both pathway table and gene table and creates one list including all pathways with all and unique genes
pw_list <- pathwaygenelist(pw_id_table_filtered, pw_gene_table,  print = FALSE)

#removes pathways with a number of genes < a certain threshold
pw_list.copy.filt <- removesmallpwinlist(pw_list, mng, print=FALSE)

#------------------Prepare for base clustering----------------

#in case either list or dataframe was created
if(exists("pw_dataframe.copy.filt")){
  clustering_height = dim(pw_dataframe.copy.filt)[1]
}else if(exists("pw_list.copy.filt")){
  clustering_height = length(pw_list.copy.filt)
}

print("preprocessing finished")

#------------------Determine k----------------

kmatrix_indices_euclidean = matrix(data = NA,nrow = clustering_height,ncol = length(indices_clustering_validation)+1)
rownames(kmatrix_indices_euclidean) <- names(pw_list.copy.filt)
colnames(kmatrix_indices_euclidean) <- c(indices_clustering_validation, "silhouette")

kmatrix_indices_spearman = matrix(data = NA,nrow = clustering_height,ncol = length(indices_clustering_validation)+1)
rownames(kmatrix_indices_spearman) <- names(pw_list.copy.filt)
colnames(kmatrix_indices_spearman) <- c(indices_clustering_validation, "silhouette")

kmatrix_indices_pearson = matrix(data = NA,nrow = clustering_height,ncol = length(indices_clustering_validation)+1)
rownames(kmatrix_indices_pearson) <- names(pw_list.copy.filt)
colnames(kmatrix_indices_pearson) <- c(indices_clustering_validation, "silhouette")


for(i in 1:clustering_height){
  
  #extract expression data for genes/pathways
  pathwayid = names(pw_list.copy.filt[i])
  print(pathwayid)
  exprdataframe <- extractexprdata(exprdata, pw_list.copy.filt, input="list", samples="all", pathwayid, FALSE)
  
  #remove NA's
  exprdataframe <- exprdataframe[complete.cases(exprdataframe),]
  exprdataframe.transposed <- t(exprdataframe)
  
  save(exprdataframe.transposed, file= paste(pathrda, subDir_gep, dataset, "_exprdataframe_transposed_",i,"_minpw",mng,".rda", sep =""))

  #compute indices for best k and save indices in matrix
  kmatrix_indices_euclidean[i,] <- indexbestk(similaritymatrix=as.matrix(exprdataframe.transposed), dissimilaritymatrix = NULL, dist="euclidean", nbmin=2, nbmax, indices_clustering_validation)

  dissimilaritymatrix_spearman <- computedistancematrix(as.matrix(exprdataframe.transposed), metric="spearman")
  kmatrix_indices_spearman[i,] <- indexbestk(similaritymatrix=as.matrix(exprdataframe.transposed), dissimilaritymatrix = dissimilaritymatrix_spearman, dist=NULL, nbmin=2, nbmax, indices_clustering_validation)

  dissimilaritymatrix_pearson <- computedistancematrix(as.matrix(exprdataframe.transposed), metric="pearson")
  kmatrix_indices_pearson[i,] <- indexbestk(similaritymatrix=as.matrix(exprdataframe.transposed), dissimilaritymatrix = dissimilaritymatrix_pearson, dist=NULL, nbmin=2, nbmax, indices_clustering_validation)
  print("determined best k - clustering 1")
  
  
}
# compute MSE over columns of kmatrix
mse_euclidean <- indexvalidation(kmatrix_indices_euclidean, actual_nb_clusterings)
mse_spearman <- indexvalidation(kmatrix_indices_spearman, actual_nb_clusterings)
mse_pearson <- indexvalidation(kmatrix_indices_pearson, actual_nb_clusterings)
print("calculated mse - clustering 1")


save(kmatrix_indices_euclidean, file= paste(pathrda, subDir_k, dataset, "_kmatrix_indices_euclidean_minpw5.rda", sep =""))
save(mse_euclidean, file= paste(pathrda, subDir_k, dataset,"_mse_euclidean_minpw5.rda", sep =""))

save(kmatrix_indices_spearman, file= paste(pathrda, subDir_k, dataset, "_kmatrix_indices_spearman_minpw5.rda", sep =""))
save(mse_spearman, file= paste(pathrda, subDir_k, dataset,"_mse_spearman_minpw5.rda", sep =""))

save(kmatrix_indices_pearson, file= paste(pathrda, subDir_k, dataset, "_kmatrix_indices_pearson_minpw5.rda", sep =""))
save(mse_pearson, file= paste(pathrda, subDir_k, dataset,"_mse_pearson_minpw5.rda", sep =""))


#------------------Temporary storage----------------

save(pw_list.copy.filt, file= paste(pathrda, subDir_gep, dataset, "_pw_list_copy_filt_minpw",mng,".rda", sep =""))
save(clustering_height, file= paste(pathrda, subDir_gep, dataset, "_clustering_height_minpw",mng,".rda", sep =""))

#---------------Load Data-----------------


load(paste(pathrda, subDir_gep, dataset, "_pw_list_copy_filt_minpw",mng,".rda", sep =""))
load(paste(pathrda, subDir_gep, dataset, "_clustering_height_minpw",mng,".rda", sep =""))


#dataframe for clustering, width: number of samples, height: number of clusterings/pathways
clustermatrixs <- data.frame(matrix(ncol = dim(exprdata)[2], nrow = clustering_height))
rownames(clustermatrixs) <- names(pw_list.copy.filt)
colnames(clustermatrixs) <- colnames(exprdata)

clustermatrixp <- data.frame(matrix(ncol = dim(exprdata)[2], nrow = clustering_height))
rownames(clustermatrixp) <- names(pw_list.copy.filt)
colnames(clustermatrixp) <- colnames(exprdata)

clustermatrixe <- data.frame(matrix(ncol = dim(exprdata)[2], nrow = clustering_height))
rownames(clustermatrixe) <- names(pw_list.copy.filt)
colnames(clustermatrixe) <- colnames(exprdata)


#cluster with a given k=2...15
for(ktest in 2:nbmax){
  #for every pathway
  for(i in 1:clustering_height){
    
    load(paste(pathrda, subDir_gep, dataset, "_exprdataframe_transposed_",i,"_minpw",mng,".rda", sep =""))
    dissimilaritymatrix_spearman <- computedistancematrix(as.matrix(exprdataframe.transposed), metric="spearman")
    dissimilaritymatrix_pearson <- computedistancematrix(as.matrix(exprdataframe.transposed), metric="pearson")
    dissimilaritymatrix_euclid <- computedistancematrix(as.matrix(exprdataframe.transposed), metric="euclid")
    
    #first clustering
    clustermatrixs[i,] <- pam(dissimilaritymatrix_spearman, ktest, diss = TRUE, metric = NULL, medoids = NULL, stand = TRUE)$clustering
    clustermatrixp[i,] <- pam(dissimilaritymatrix_pearson, ktest, diss = TRUE, metric = NULL, medoids = NULL, stand = TRUE)$clustering
    clustermatrixe[i,] <- pam(dissimilaritymatrix_euclid, ktest, diss = TRUE, metric = NULL, medoids = NULL, stand = TRUE)$clustering
    
  }
  print("first clustering complete")
  
  #calculate voting and similarity matrix
  similaritymatrixs <- samplesimilarity(clustermatrixs)
  similaritymatrixp <- samplesimilarity(clustermatrixp)
  similaritymatrixe <- samplesimilarity(clustermatrixe)
  print(similaritymatrixe)
  
  print("Calculated Similarity matrix.")
  
  save(clustermatrixs, file= paste(pathrda, subDir_cc ,dataset,"_clustermatrixs_k",ktest,"_minpw",mng,".rda", sep =""))
  save(clustermatrixp, file= paste(pathrda, subDir_cc ,dataset,"_clustermatrixp_k",ktest,"_minpw",mng,".rda", sep =""))
  save(clustermatrixe, file= paste(pathrda, subDir_cc ,dataset,"_clustermatrixe_k",ktest,"_minpw",mng,".rda", sep =""))
  
  save(similaritymatrixs, file= paste(pathrda, subDir_cc ,dataset,"_votingmatrixs_k",ktest,"_minpw",mng,".rda", sep =""))
  save(similaritymatrixp, file= paste(pathrda, subDir_cc ,dataset,"_votingmatrixp_k",ktest,"_minpw",mng,".rda", sep =""))
  save(similaritymatrixe, file= paste(pathrda, subDir_cc ,dataset,"_votingmatrixe_k",ktest,"_minpw",mng,".rda", sep =""))
  
  print("Saved matrix.")
  
  #-------------spearman------------
  dissimilaritymatrix_spearman <- as.dist(1-similaritymatrixs/clustering_height)
  
  print("calculated distance matrix 2")
  save(dissimilaritymatrix_spearman, file= paste(pathrda, subDir_cc ,dataset,"_clustering1_dissimilaritymatrix_spearman_k",ktest,"_minpw",mng,".rda", sep =""))
  consensusclustering_spearman <- pam(dissimilaritymatrix_spearman, ktest, diss = TRUE)
  save(consensusclustering_spearman, file=paste(pathrda, subDir_cc ,dataset, "_consensusclustering_spearman_k",ktest,"_minpw",mng,".rda", sep =""))
  
  #-------------pearson------------
  dissimilaritymatrix_pearson <- as.dist(1-similaritymatrixp/clustering_height)
  print("calculated distance matrix 2")
  save(dissimilaritymatrix_pearson, file= paste(pathrda, subDir_cc ,dataset,"_clustering1_dissimilaritymatrix_pearson_k",ktest,"_minpw",mng,".rda", sep =""))
  consensusclustering_pearson <- pam(dissimilaritymatrix_pearson, ktest, diss = TRUE)
  save(consensusclustering_pearson, file=paste(pathrda, subDir_cc ,dataset, "_consensusclustering_pearson_k",ktest,"_minpw",mng,".rda", sep =""))
  
  #-------------euclid------------
  dissimilaritymatrix_euclid <- as.dist(1-similaritymatrixe/clustering_height)
  print("calculated distance matrix 2")
  save(dissimilaritymatrix_euclid, file= paste(pathrda, subDir_cc ,dataset,"_clustering1_dissimilaritymatrix_euclid_k",ktest,"_minpw",mng,".rda", sep =""))
  consensusclustering_euclid <- pam(dissimilaritymatrix_euclid, ktest, diss = TRUE)
  save(consensusclustering_euclid, file=paste(pathrda, subDir_cc ,dataset, "_consensusclustering_euclid_k",ktest,"_minpw",mng,".rda", sep =""))
  
}
print("Consensus Clustering is finished.")