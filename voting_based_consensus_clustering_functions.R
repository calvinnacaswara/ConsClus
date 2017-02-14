
# Freie Universitaet Berlin
# Fachbereich Mathematik und Informatik
# Institut fuer Informatik, Institut fuer Bioinformatik

# Calvinna Yosephine Caswara

#------------------Packages--------------- 
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
library(flexclust)
library(corrplot)

#------------------Functions------------------

#############################################

# FUNCTION rowvar
# computes variance over rows of a matrix

#############################################

rowvar <- function(x,percentage) {
  a = rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
  xnew = tail(sort(a),round(length(a)*(percentage/100)))
  return(names(xnew))
}


#############################################

# FUNCTION excludepw
# excludes given pathways that are too big or irrelevant

#############################################

excludepw <- function(pw_id_table = NULL,
                      elist = NULL,
                      print = FALSE){
  
  
  #convert table into vector
  pw_id_table_filt <- unlist(pw_id_table[,1, drop=FALSE]) 
  
  #go through pathway list
  for (i in length(pw_id_table_filt):1){
    #go through list with to be excluded pathways
    for (j in length(elist):1){
      
      #comparison (if not null)
      if (!is.na(as.character(pw_id_table_filt[i]))){
        if (as.character(pw_id_table_filt[i]) == as.character(elist[j])){
          
          if (print == TRUE){
            #cat("removing: ", pw_id_table_filt[i])
            #print(pw_id_table_filt[i])
            #cat("length of table: ", length(pw_id_table_filt))
          }
          #removing element
          pw_id_table_filt <- pw_id_table_filt[-i]
        }
      }
    }
  }
  print("Removed chosen pathways from list.")
  return(pw_id_table_filt)
}



#############################################

# FUNCTION pathwaygenelist
# reads both pathway table and gene table and creates one list 
# including all pathways with all and unique genes

#############################################

pathwaygenelist <- function(pw_id_table = NULL,
                            pw_gene_table = NULL,
                            print = FALSE){
  
  #initializing
  index_gene_match <- 1
  list_of_pw = vector("list", length(pw_id_table))
  table_pw_genes = vector("list", length(pw_id_table))
  pwlist <- list()
  
  
  pwnames = unique(gsub(":.*$", "", pw_id_table))
  
  #go through pathway list
  for (pw in pwnames) {
    
    pwsearchstring = paste("^",pw,":",sep="")
    
    genestring = pw_gene_table[grep(pwsearchstring, pw_gene_table[,1]),2]
    
    #split gene string after ~
    #print("Splitting gene string")
    gene_list.split <- c(strsplit(genestring, "~"))
    
    #make 1 vector of all vectors
    gene_list.split.unl <- c(unlist(gene_list.split))
    
    #remove empty strings
    #print("Filtering gene list")
    gene_list.split.unl.filt <- gene_list.split.unl[gene_list.split.unl != ""]
    
    #make sure there is no NA
    lapply(gene_list.split.unl.filt, function(x) x[!is.na(x)])
    
    #only keep uniques
    #print("Keeping only unique genes")
    gene_list.split.unl.filt.uni <- unique(gene_list.split.unl.filt)
    
    if(print==TRUE){
      cat("Pathway: ", pw, "\n")
      print(gene_list.split.unl.filt.uni)
    }
    
    pwlist[[pw]] <- gene_list.split.unl.filt.uni
    
  }
  return(pwlist)
}

#############################################
# ------------------------------------
# NEW FUNCTION pathwaygenelist2
# reads both pathway table and gene table and creates one list 
# including all pathways with all and unique genes

#############################################

pathwaygenelist2 <- function(pw_id_table = NULL,
                             pw_gene_table = NULL, 
                             df_length = 200,
                             print = FALSE){
  
  #initializing
  index_gene_match <- 1
  list_of_pw = c()
  table_pw_genes = c()
  process_df = ""
  pw_dataframe = data.frame(matrix(ncol = df_length, nrow = 1))
  pw_datalist = c()
  
  
  #go through pathway list
  for (i in 1:length(pw_id_table)){
    
    collect_sec = c()
    
    #go through gene list
    for (j in index_gene_match:dim(pw_gene_table)[1]){
      
      #apply strsplit on every row (2nd argument: 1) and split after ":"
      genes_id <- sapply(strsplit(c(pw_gene_table[j,1]), ":"),"[", 1)
      
      #comparison with the filtered pathway list
      if (pw_id_table[i] == genes_id){
        
        #initialize new
        collect_first <- c()
        
        #split gene string after ~
        gene_list <- c(sapply(strsplit(c(pw_gene_table[j,2]), "~"),"["))
        
        #remove empty strings "" and NA
        collect_first <- c(unlist(gene_list[gene_list != ""]))
        lapply(collect_first, function(x) x[!is.na(x)])
        
        collect_sec <- c(collect_sec, collect_first)
        
        #save index of last gene matched with pathway id to save time
        index_gene_match <- j+1
        
        #set to its pathway id to remind there is a vector to be processed
        process_df <- "pw_id_table[i]"
      }
      
    } 
    if (process_df != ""){
      
      #only keep unique genes
      collect_row <- unique(unlist(collect_sec))
      collect_row <- unique(collect_row)
      
      #convert into dataframe
      a <- length(collect_row)
      while(a < df_length){
        collect_row <- c(collect_row, "empty")
        a = a+1
      }
      
      pw_dataframe <- rbind(pw_dataframe, collect_row)
      #print(dim(pw_dataframe))
      process_df = ""
      collect_row <- c()
    }
    
  }
  
  #now delete first row (because NA)
  pw_dataframe <- pw_dataframe[-1,]
  
  #set pathway names
  rownames(pw_dataframe) <- pw_id_table
  
  return(pw_dataframe)
}



#############################################

# FUNCTION removesmallpwindataframe
# removes pathways with a number of genes < a certain threshold

#############################################

removesmallpwindataframe <- function(list_of_pw = NULL, minimal_gene_nb_threshold = 10, print = FALSE){
  
  #counter for removed pathways
  removed_counter <- 0
  #empty list
  index_remove <- 0
  
  for (i in 1:dim(list_of_pw)[1]){
    if(as.integer(length(list_of_pw[i,])-(sum(list_of_pw[i,] == ""))) < as.integer(minimal_gene_nb_threshold)){
      
      #save index of pathway (to be removed) in list
      index_remove <- c(index_remove, i)
      removed_counter = removed_counter+1
      
    }
  }
  
  #deletion 
  list_of_pw <- list_of_pw[-(index_remove),]
  
  #print(dim(list_of_pw))
  
  #print("-----------------------------------")
  #cat("removed pathways with <",minimal_gene_nb_threshold, " genes: ", removed_counter)
  
  return(list_of_pw)
}

#############################################

# FUNCTION removesmallpwindatalist
# removes pathways with a number of genes < a certain threshold

#############################################

removesmallpwinlist <- function(list_of_pw = NULL, minimal_gene_nb_threshold = 10, print = FALSE){
  
  #counter for removed pathways
  removed_counter <- 0
  #empty list
  index_remove <- 0
  
  for (i in 1:length(list_of_pw)){
    if(length(list_of_pw[[i]]) < minimal_gene_nb_threshold){
      
      #save index of pathway (to be removed) in list
      index_remove <- c(index_remove, i)
      removed_counter = removed_counter+1
      
    }
  }
  
  #cat("remove elements with index: ", index_remove, "\n")
  #deletion 
  list_of_pw[index_remove] <- NULL
  
  #cat("final length of list: ", length(list_of_pw), "\n")
  
  print("-----------------------------------")
  cat("removed pathways with <",minimal_gene_nb_threshold, " genes: ", removed_counter)
  
  return(list_of_pw)
}

#############################################

# FUNCTION extractexprdata
# extracts expression data of every gene (in certain/all samples) in a certain pathway
# works both for list and dataframe as input 

#############################################

extractexprdata <- function(exprtable = NULL, pw_list = NULL, input="", samples="all", pathwayid = NULL, print = FALSE){
  
  
  if(input!="list" && input!="dataframe"){
    print("Please specify input type of pathway_list: list, dataframe")
  }
  
  #input type list
  else if(input=="list"){
    
    #unlist list, select rows in expression dataframe with geneids
    #geneids are saved in pw_list (with their pathwayid as name)
    if(samples=="all"){
      
      #select genes in dataframe (pathwayid, geneid)
      pathway_gene_row = pw_list[[pathwayid]]
      
      #do not select "empty" fields
      pathway_gene_row_noempt <- pathway_gene_row[pathway_gene_row !=""]
      
      #check if gene exists in dataframe
      #for(elem in pathway_gene_row_noempt){
      #  if(!(elem %in% rownames(exprtable))){
      #    pathway_gene_row_noempt_exist <- pathway_gene_row_noempt[pathway_gene_row_noempt != elem]
      #  }
      #}
      #print(pathway_gene_row_noempt_exist)
      
      
      inds <- which(rownames(exprtable) %in% pathway_gene_row_noempt)
      
      #select expression data in dataframe (geneid, exprdata per sample)
      exprdataframe <- exprtable[inds,]
      
    }
    else{
      exprdataframe <- exprtable[pathway_gene_row_noempt, samples]
    }
    
    return(exprdataframe)
    
  }
  
  #input type dataframe
  else{
    
    #unlist list, select rows in expression dataframe with geneids
    #geneids are saved in pw_list (with their pathwayid as name)
    if(samples=="all"){
      
      #select genes in dataframe (pathwayid, geneid)
      pathway_gene_row = pw_list[pathwayid,]
      
      #do not select "empty" fields
      pathway_gene_row_noempt <- pathway_gene_row[pathway_gene_row !=""]
      
      inds <- which(rownames(exprtable) %in% pathway_gene_row_noempt)
      
      
      #select expression data in dataframe (geneid, exprdata per sample)
      exprdataframe <- exprtable[inds,]
      
    }
    else{
      exprdataframe <- exprtable[pathway_gene_row_noempt, samples]
    }
    
    return(exprdataframe)
  }
  
  
}


#############################################

# FUNCTION list2dataframe
# makes a dataframe out of a list with minimal size

#############################################

list2dataframe <- function(pwlist = NULL,
                           print = FALSE){
  max.nb.genes = 0
  
  for(pw in pwlist){
    if(length(pw)>max.nb.genes){
      max.nb.genes = length(pw)
    }
  }
  
  pw_dataframe = data.frame(matrix(ncol = max.nb.genes, nrow = 1))
  pw_datalist = c()
  
  
  for(pw in pwlist){
    #convert into dataframe
    a <- length(pw)
    while(a < max.nb.genes){
      pw <- c(pw, "")
      pw_datalist <- pw
      a = a+1
    }
    pw_dataframe <- rbind(pw_dataframe, pw_datalist)
  }
  
  #now delete first row (because NA)
  pw_dataframe <- pw_dataframe[-1,]
  
  #set pathway names
  rownames(pw_dataframe) <- names(pwlist)
  
  #print(dim(pw_dataframe))
  return(pw_dataframe)
}



#############################################

# FUNCTION samplesimilarity
# counts how many times a sample pair was clustered in 
# single clusters and divides it through the number of 
# clusterings or pathways

#############################################


samplesimilarity <- function(clusterings = NULL, print = FALSE){
  
  
  #create n*n matrix (pairwise voting)
  votingmatrix <- data.frame(matrix(ncol = dim(clusterings)[2], nrow = dim(clusterings)[2]))
  rownames(votingmatrix) <- colnames(clusterings)
  colnames(votingmatrix) <- colnames(clusterings)
  
  
  #fill matrix with zeros
  votingmatrix[is.na(votingmatrix)] <- 0
  
  #distance between first voting sample and second sample
  k = 1
  
  #for every row/clustering/pathway
  
  loop <- combn(c(1:(dim(clusterings)[2])),2)
  
  for(i in 1:dim(clusterings)[1]){
    
    #index of second sample
    for (h in 1:(dim(loop)[2])){
      k = loop[2,h]
      j = loop[1,h]
      
      if(clusterings[i,j]==clusterings[i,k]){
        votingmatrix[j,k] <- votingmatrix[j,k] + 1 
        votingmatrix[k,j] <- votingmatrix[k,j] + 1
      }
    }
    print(i)
  }
  
  #create similarity matrix by dividing every observation/voting through the number of clusterings/pathways
  #votingmatrix <- votingmatrix/dim(clusterings)[1]
  print(votingmatrix)
  
  return(votingmatrix)
}


#############################################

# FUNCTION indexbestk
# computes different (non graphical) indices for best k

#############################################


indexbestk <- function(similaritymatrix = NULL, dissimilaritymatrix = NULL, dist = NULL, nbmin = 2, nbmax = 15, indices_clustering_validation = NULL, print = FALSE){
  
  indices <- c()  
  
  for(i in 1:length(indices_clustering_validation)){
    indices[i] <- NbClust(similaritymatrix, diss=dissimilaritymatrix, distance=dist, min.nc=nbmin, max.nc=nbmax,
                          method="centroid", index=(indices_clustering_validation[i]))$Best.nc[1]
    
    #zum testen
    #indices[i] <- NbClust(as.matrix(exprdataframe.transposed), diss = NULL, distance="euclidean", min.nc=2, max.nc=23, method="centroid",indices_clustering_validation[i])$Best.nc[1]
  }
  indices[length(indices_clustering_validation)+1] <- pamk(similaritymatrix, krange=2:15, criterion="asw", usepam = TRUE)$nc
  #print(indices)
  
  return(indices)
  
}

#############################################

# FUNCTION indexvalidation
# computes mean square error (MSE) for every column/sample

#############################################


indexvalidation <- function(kmatrix_indices = NULL, actual_nb_clusterings = NULL, print = FALSE){
  
  mse_index <- list()  
  
  #go through columns
  for(j in 1:dim(kmatrix_indices)[2]){
    #calculate mse for every column/sample
    mse_index[j] <- (1/dim(kmatrix_indices)[1])*(sum((kmatrix_indices[,j]-actual_nb_clusterings)^2))
  }
  names(mse_index) <- colnames(kmatrix_indices)
  print(mse_index)
  #sort mse of indices and select the first 5 elements (take the top 5 indices with the smalles mse)
  mse_sorted = sort(unlist(mse_index))
  
  
  return(mse_sorted)
  
}

#############################################

# FUNCTION computedistancematrix
# computes distance matrix 
# input: data (as matrix)

#############################################

computedistancematrix <- function(data = NULL, metric = "none", print = FALSE){
  
  #compute dissimilarity matrix
  if (metric == "pearson"){
    dissimilaritymatrix <- cor.dist(as.matrix(data))
  }else if(metric=="spearman"){
    dissimilaritymatrix <- spearman.dist(as.matrix(data))
  }else if(metric=="euclid"){
    dissimilaritymatrix <- euc(as.matrix(data))
  }else{
    dissimilaritymatrix <- 1-data
    print("Neither spearman nor pearson correlation distance were chosen, so the normal dissmilarity matrix was computed.")
  }
  return (dissimilaritymatrix)
}

#############################################

# FUNCTION sdcheck
#check if standard deviation > 0

#############################################

sdcheck <- function(similaritymatrix = NULL, print = FALSE){
  null_sd = FALSE
  row_standard_deviation = rowSds(as.matrix(similaritymatrix), na.rm=TRUE)
  for (elem in row_standard_deviation){
    if (elem == 0){
      null_sd = TRUE
    }
  }
  return(null_sd)
}

#############################################

# FUNCTION randindex
# calculates measure of the similarity between two data clusterings

#############################################

randindex <- function(clustering1 = NULL, clustering2 = NULL,print = FALSE){
  
  a=0
  b=0
  c=0
  d=0
  
  loop <- combn(c(1:length(clustering1)),2)
  for (h in 1:(dim(loop)[2])){
    j = loop[1,h]
    k = loop[2,h]
    
    if(clustering1[j]==clustering1[k] && clustering2[j]==clustering2[k]){
      a = a+1
    }
    else if(clustering1[j]==clustering1[k] && clustering2[j]!=clustering2[k]){
      d = d+1
    }
    else if(clustering1[j]!=clustering1[k] && clustering2[j]==clustering2[k]){
      c = c+1
    }
    else{
      b = b+1
    }
  }
  
  RI = (a+b)/(a+b+c+d)
  return(RI)
}


#############################################

# FUNCTION valmatagg

#############################################

valmatagg <- function(path, dataset = NULL, distance = NULL,consensusclustering_original, mng, print = FALSE){
  
  #create matrix
  valmatrix <- data.frame(matrix(ncol = 14, nrow = 3))
  rownames(valmatrix) <- c("RI","ARI","NMI")
  colnames(valmatrix) <- c(2:15)
  
  for(i in 1:14){
    if(distance=="euclid"){
      load(paste(path, dataset,"_consensusclustering_euclid_k",i+1,"_minpw",mng,".rda", sep=""))
      cons = consensusclustering_euclid$clustering
    }else if(distance=="spearman"){
      load(paste(path, dataset,"_consensusclustering_spearman_k",i+1,"_minpw",mng,".rda", sep=""))
      cons = consensusclustering_spearman$clustering
    }else if(distance=="pearson"){
      load(paste(path, dataset,"_consensusclustering_pearson_k",i+1,"_minpw",mng,".rda", sep=""))
      cons = consensusclustering_pearson$clustering
    }
    
    
    valmatrix[1,i] = randindex(cons, consensusclustering_original)
    valmatrix[2,i] = randIndex(cons, consensusclustering_original)
    valmatrix[3,i] = cl_agreement(as.cl_ensemble(list(as.cl_partition(cons),as.cl_partition(consensusclustering_original))),method="NMI")
    
  }
  return(valmatrix)
}
