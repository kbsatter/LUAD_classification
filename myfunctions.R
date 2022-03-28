
# Analysis Date: June 7, 2019
# Author: Paul Tran
# Title: Function to test multiple UMAP paramters




see<-function(df, row = 10, column = 10){
  df[1:row,1:column]
}


################ DBU #################
umap_parameters_gridsearch <-
  function(dataset, genelist, colorlabel, metric = "euclidean") {
    "
Searches number of neighbors 2-50, minimum distance 0-0.99, and number of random genes 10-10000
Outputs pdf with 
"
    require(umap)
    require(ggplot2)
    custom.config = umap.defaults
custom.config$metric = metric
count =0


for (i in c(5, 20, 40)){
  for (j in c(0.001, 0.01, 0.1, 0.25)){
    for (k in c(1000,2500, 5000)){
      # i=5
      # j=0.0
      # k=250
      print(paste("Metric=",metric,"N_neighor = ",i,"Min_dist=",j,"Gene_no",k))
      custom.config$n_neighbors = i
    custom.config$min_dist = j
    rand_X_genes<-sample(genelist,k)
    brain.umap.2 = umap(dataset[,rand_X_genes], custom.config)
    count = count + 1
    if (count == 1){
      g<-ggplot(data.frame(brain.umap.2$layout),aes(X1,X2, color=colorlabel))+
        geom_point() +
        theme(legend.position = c(0.9, 0.5)) +
        ggtitle(paste("n_neighbors=",i, "min_dist=",j, "no_genes=",k))
      plot(g)
    }
    else{    
      g<-ggplot(data.frame(brain.umap.2$layout),aes(X1,X2, color=colorlabel))+
        geom_point() +
        theme(legend.position= "none")+
        ggtitle(paste("n_neighbors=",i, "min_dist=",j, "no_genes=",k))
      plot(g)
    }
    }
  }
}

}


umapiterations<-
  function(dataset, 
           genes,
           no_iters = 10,
           rand_gene_no = 5000, 
           n_neighbors = 20, 
           min_dist = 0.1, 
           metric = "manhattan",
           samp_prop = 0.5
           )
  {
    "parallelized UMAP iterations"
    require(foreach)
    
    #
    if(rand_gene_no>length(genes)){paste("Too many random genes for number of genes given")}
    
    ##parameters
    custom.config = umap::umap.defaults
    no_iters=no_iters
    rand_gene_no=rand_gene_no
    custom.config$n_neighbors = n_neighbors
    custom.config$min_dist = min_dist
    custom.config$metric = metric
    
    # samp_matrix<-caret::createDataPartition(study2,p=samp_prop,times = no_iters)
    # samp_matrix<-do.call(cbind.data.frame, samp_matrix)
    # umap_position_matrix<-bigstatsr::FBM(dim(samp_matrix)[1],2*no_iters)
    
    umap_position_matrix<-bigstatsr::FBM(dim(dataset)[1],2*no_iters)
    
    start<-Sys.time()
    
    cl <- parallel::makeForkCluster(bigstatsr::nb_cores())
    doParallel::registerDoParallel(cl)

    foreach (i = 1:no_iters, .combine = "c") %dopar% {
      # my.umap = umap::umap(dataset[samp_matrix[,i],sample(genes,rand_gene_no)],custom.config)
      set.seed(252+i*3)
      my.umap = umap::umap(dataset[,sample(genes,rand_gene_no)],custom.config)
      
      x<-(i*2)-1
      umap_position_matrix[,x]<-my.umap$layout[,1]
      umap_position_matrix[,x+1]<-my.umap$layout[,2]
      
      NULL
    }
    parallel::stopCluster(cl)
    
    stop<-Sys.time()
    print(stop-start)
    
    # return(list(umap_position_matrix,samp_matrix))
    umap_position_matrix
  }


umapiterations2<-
  function(dataset, 
           genes,
           no_iters = 10,
           rand_gene_no = 5000, 
           n_neighbors = 20, 
           min_dist = 0.1, 
           metric = "manhattan")
  {
    "non-parallelized umap iterations"
    ##parameters##
    custom.config = umap::umap.defaults
    no_iters=no_iters
    rand_gene_no=rand_gene_no
    custom.config$n_neighbors = n_neighbors
    custom.config$min_dist = min_dist
    custom.config$metric = metric
    
    umap_position_matrix<-matrix(NA,dim(dataset)[1],2*no_iters)
    start<-Sys.time()

    for (i in 1:no_iters) {
      # i=1
      print(paste("Iteration #",i))
      my.umap = umap::umap(dataset[,sample(genes,rand_gene_no)],custom.config)
      x<-(i*2)-1
      umap_position_matrix[,x]<-my.umap$layout[,1]
      umap_position_matrix[,x+1]<-my.umap$layout[,2]
      
      NULL
    }
    
    stop<-Sys.time()
    print(stop-start)
    umap_position_matrix
    
  }



dbutest<-
  function(umap_position_matrix,iter,k = 15,eps = 0.5){
    "explore optimal parameters for dbscan algorithm using data from UMAP iterations"
    # par(mfrow = c(2, 1))
    x<-(iter*2)-1
    # Compute DBSCAN using fpc package
    df<-umap_position_matrix[,x:(x+1)]
    #Method for determining the optimal eps value
    a<-dbscan::kNNdistplot(df, k = k)
    abline(h = eps, lty = 2)
    db <- fpc::dbscan(df, eps = eps, MinPts = k)
    b<-plot(db, df, frame = FALSE,xlab="",ylab="")
    return(list(a, b))
  }


dbuiters<-
  function(umap_position_matrix,k=15, eps=0.5, filename="DBU_iters_"){
    "run dbscan on all UMAP iterations based on optimezed dbscan parameters from dbutest function"
    no_iters<-dim(umap_position_matrix)[2]/2
    
    pdf(paste(filename, "Reps",no_iters,".pdf"))
    par(mfrow=c(4,3),mar=c(0,0,0,0))
    db_umap_res<-matrix(NA,nrow = dim(umap_position_matrix)[1],ncol = no_iters)

        for (i in 1:no_iters){
          x<-(i*2)-1
          df<-umap_position_matrix[,x:(x+1)]
          db <- fpc::dbscan(df, eps = eps, MinPts = k)
      
          plot(db, df, frame = FALSE,xlab="",ylab="")
          legend('topleft',paste(i),bty = 'n',text.col='red',text.font = 2)
      
          db_umap_res[,i]<-db$cluster
          print(i)
    }
    dev.off()
    
    colnames(db_umap_res)<-1:no_iters
    
    db_umap_res
  }

change_grp_names5<-function(x){
  require(dplyr)
  "takes single DBU iteration and renumbers the most common group as 1,
  second most common as 2, and so on for the top 4 groups, all others
  are 0"
  x[x==names(sort(table(x)))[length(unique(x))]]<-"Grp1"
  x[x==names(sort(table(x)))[length(unique(x))-1]]<-"Grp2"
  x[x==names(sort(table(x)))[length(unique(x))-2]]<-"Grp3"
  x[x==names(sort(table(x)))[length(unique(x))-3]]<-"Grp4"
  x[x==names(sort(table(x)))[length(unique(x))-4]]<-"Grp5"
  #x[x!="Grp1"|x!="Grp2"|x!="Grp3"|x!="Grp4"]<-NA
  x<-recode(x, Grp1 = 1,
            Grp2 = 2,
            Grp3 = 3,
            Grp4 = 4,
            Grp5 = 5)
}


change_grp_names4<-function(x){
  require(dplyr)
  "takes single DBU iteration and renumbers the most common group as 1,
  second most common as 2, and so on for the top 4 groups, all others
  are 0"
  x[x==names(sort(table(x)))[length(unique(x))]]<-"Grp1"
  x[x==names(sort(table(x)))[length(unique(x))-1]]<-"Grp2"
  x[x==names(sort(table(x)))[length(unique(x))-2]]<-"Grp3"
  x[x==names(sort(table(x)))[length(unique(x))-3]]<-"Grp4"
  #x[x!="Grp1"|x!="Grp2"|x!="Grp3"|x!="Grp4"]<-NA
  x<-recode(x, Grp1 = 1,
            Grp2 = 2,
            Grp3 = 3,
            Grp4 = 4)
}

change_grp_names3<-function(x){
  require(dplyr)
  "takes single DBU iteration and renumbers the most common group as 1,
  second most common as 2, and so on for the top 3 groups, all others
  are 0"
  x[x==names(sort(table(x)))[length(unique(x))]]<-"Grp1"
  x[x==names(sort(table(x)))[length(unique(x))-1]]<-"Grp2"
  x[x==names(sort(table(x)))[length(unique(x))-2]]<-"Grp3"
  x[x!="Grp1"|x!="Grp2"|x!="Grp3"]<=0
  x<-recode(x, Grp1 = 1,
            Grp2 = 2,
            Grp3 = 3)
}


makegrpcalls<-
  function(dbu_clean, cutoff = 75){
    "takes input dbu grouping matrix from dbuiters function and outputs plurality voted groups;colna
    CLEAN the matrix first!"
#find ambiguous samples
all_maxgrp<-c()

for (j in 1:dim(dbu_clean)[1]){
  #j=4
  datab<-table(as.character(dbu_clean[j,]))
  maxgrp<-names(datab)[(which.max(datab))]
  maxperc<-max(table(as.character(dbu_clean[j,]))*100/dim(dbu_clean)[2])
  toadd<-c(maxgrp,maxperc)
  all_maxgrp<-rbind(all_maxgrp,toadd)
}
rownames(all_maxgrp)<-rownames(dbu_clean)
colnames(all_maxgrp)<-c("Group","Proportion")

ambiloc<-which(all_maxgrp[,1]=="0"|as.numeric(all_maxgrp[,2])<cutoff)
all_maxgrp[ambiloc,]->ambig

#define DBUgrp to classify all samples
DBUgrp<-all_maxgrp[,1]
DBUgrp[ambiloc]<-"Ambi"
table(DBUgrp)

all_maxgrp<-cbind.data.frame(all_maxgrp,DBUgrp)
all_maxgrp
  }



makegrpcalls_clust<-function(db_umap_res3_4,k){
  samp_clust<-hclust(dist(db_umap_res3_4,method = "manhattan"))
  samp_clust6<-cutree(samp_clust,k=k)
  
  iter_clust<-hclust(dist(t(db_umap_res3_4),method = "manhattan"))
  
  sample_cluster_diversity <- function(vec1){
    vec1<-as.vector(vec1[!is.na(vec1)])
    for (i in length(vec1):2){
      x <- abs(vec1[i]-vec1[i-1])
      if (i == length(vec1)){
        y = x
      }
      else {
        y <- y + x
      }
    }
    y
  }
  
  
  samp_clus_div_res<-apply(db_umap_res3_4[,iter_clust$order],1,sample_cluster_diversity)
  consens_calls<-cbind.data.frame(samp_clust6,samp_clus_div_res)
  consens_calls
}

connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  "
  from line 693 of M3C.R code
  https://github.com/crj32/M3C/blob/master/R/M3C.R
  "
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey 
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId
  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}

group_consensus_freq<-function(umap_position_matrix2,db_umap_res){
n=nrow(umap_position_matrix2)
mCount = matrix(c(0),ncol=n,nrow=n)

for (i in 1:ncol(db_umap_res)){
  mCount<-connectivityMatrix(db_umap_res[,i],mCount,1:n)
  print(i)
}

Mfreq<-mCount/max(mCount)
}

################ ROC/PR #####################
multi_ROC_5grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G5_true=fastDummies::dummy_cols(Lim_obs)[,6],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4,
    G5_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP5)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)

  roc_test$AUC
}


multi_ROC_4grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)
  
  roc_test$AUC
}


multi_ROC_3grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)
  
  roc_test$AUC
}

plot_multi_ROC_5grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G5_true=fastDummies::dummy_cols(Lim_obs)[,6],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4,
    G5_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP5)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)
  plot_roc_df<-multiROC::plot_roc_data(roc_test)
  torm<-which(plot_roc_df$Group=="Micro"|plot_roc_df$Group=="Macro")
  plot_roc_df$Group<-dplyr::recode(plot_roc_df$Group, 
                                   G1 = "LUAD1", 
                                   G2 = "LUAD2", 
                                   G3 = "LUAD3", 
                                   G4 = "LUAD4", 
                                   G5 = "LUAD5")
  
  g<-ggplot(plot_roc_df[-torm,], aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF",
                                "#5DC863FF",
                                "#FDE725FF"))+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  # ggsave("Figure4REMBRANDT_ROC.png",
  #        width = 3,
  #        height = 2.5,
  #        dpi = 300,
  #        units = "in")
  g
  
}


plot_multi_ROC_4grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)
  plot_roc_df<-multiROC::plot_roc_data(roc_test)
  torm<-which(plot_roc_df$Group=="Micro"|plot_roc_df$Group=="Macro")
  plot_roc_df$Group<-dplyr::recode(plot_roc_df$Group, 
                                   G1 = "LUAD1", 
                                   G2 = "LUAD2", 
                                   G3 = "LUAD3", 
                                   G4 = "LUAD4")
  
  g<-ggplot(plot_roc_df[-torm,], aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF",
                                "#5DC863FF"))+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  # ggsave("Figure4REMBRANDT_ROC.png",
  #        width = 3,
  #        height = 2.5,
  #        dpi = 300,
  #        units = "in")
  g
  
}


plot_multi_ROC_3grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3)
  
  
  
  roc_test <- multiROC::multi_roc(roc_data)
  plot_roc_df<-multiROC::plot_roc_data(roc_test)
  torm<-which(plot_roc_df$Group=="Micro"|plot_roc_df$Group=="Macro")
  plot_roc_df$Group<-dplyr::recode(plot_roc_df$Group, 
                                   G1 = "LUAD1", 
                                   G2 = "LUAD2", 
                                   G3 = "LUAD3")
  
  g<-ggplot(plot_roc_df[-torm,], aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF"))+
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  # ggsave("Figure4REMBRANDT_ROC.png",
  #        width = 3,
  #        height = 2.5,
  #        dpi = 300,
  #        units = "in")
  g
  
}

multi_PR_5grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G5_true=fastDummies::dummy_cols(Lim_obs)[,6],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4,
    G5_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP5)
  
  
  pr_res <- multiROC::multi_pr(roc_data)
  pr_res$AUC
}


multi_PR_4grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4)
  
  
  pr_res <- multiROC::multi_pr(roc_data)
  pr_res$AUC
}


multi_PR_3grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3)
  
  
  pr_res <- multiROC::multi_pr(roc_data)
  pr_res$AUC
}


plot_multi_PR_5grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G5_true=fastDummies::dummy_cols(Lim_obs)[,6],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4,
    G5_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP5)
  
  pr_res <- multiROC::multi_pr(roc_data)
  plot_pr_df<-multiROC::plot_pr_data(pr_res)
  torm<-which(plot_pr_df$Group=="Micro"|plot_pr_df$Group=="Macro")
  plot_pr_df$Group<-dplyr::recode(plot_pr_df$Group, 
                                  G1 = "LUAD1", 
                                  G2 = "LUAD2", 
                                  G3 = "LUAD3", 
                                  G4 = "LUAD4",
                                  G5 = "LUAD5")
  
  ggplot(plot_pr_df[-torm,], aes(x=Recall, y=Precision)) + 
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF",
                                "#5DC863FF",
                                "#FDE725FF"))+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  
}


plot_multi_PR_4grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G4_true=fastDummies::dummy_cols(Lim_obs)[,5],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3,
    G4_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP4)
  
  pr_res <- multiROC::multi_pr(roc_data)
  plot_pr_df<-multiROC::plot_pr_data(pr_res)
  torm<-which(plot_pr_df$Group=="Micro"|plot_pr_df$Group=="Macro")
  plot_pr_df$Group<-dplyr::recode(plot_pr_df$Group, 
                                  G1 = "LUAD1", 
                                  G2 = "LUAD2", 
                                  G3 = "LUAD3", 
                                  G4 = "LUAD4")
  
  ggplot(plot_pr_df[-torm,], aes(x=Recall, y=Precision)) + 
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF",
                                "#5DC863FF"))+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  
}


plot_multi_PR_3grp<-function(Lim_obs,model,Lim_exp){
  roc_data<-data.frame(
    G1_true=fastDummies::dummy_cols(Lim_obs)[,2],
    G2_true=fastDummies::dummy_cols(Lim_obs)[,3],
    G3_true=fastDummies::dummy_cols(Lim_obs)[,4],
    G1_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP1,
    G2_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP2,
    G3_pred_m1=predict(my_glmnet_model,Lim_exp,type = "prob")$TP3)
  
  pr_res <- multiROC::multi_pr(roc_data)
  plot_pr_df<-multiROC::plot_pr_data(pr_res)
  torm<-which(plot_pr_df$Group=="Micro"|plot_pr_df$Group=="Macro")
  plot_pr_df$Group<-dplyr::recode(plot_pr_df$Group, 
                                  G1 = "LUAD1", 
                                  G2 = "LUAD2", 
                                  G3 = "LUAD3")
  
  ggplot(plot_pr_df[-torm,], aes(x=Recall, y=Precision)) + 
    geom_path(aes(color = Group), size=1.5) +
    scale_color_manual(values=c("#440154FF", 
                                "#3B528BFF", 
                                "#21908CFF"))+
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.justification=c(1, 0), legend.position=c(.95, .05),
          legend.title=element_blank(), 
          legend.background = element_rect(fill=NULL, size=0.5, 
                                           linetype="solid", colour ="black"))
  
}

norm.interval = function(data, variance = var(data), conf.level = 0.95) {
  z = qnorm((1 - conf.level)/2, lower.tail = FALSE)
  xbar = mean(data)
  sdx = sqrt(variance/length(data))
  c(xbar - z * sdx, xbar + z * sdx)
}



fourgrpDEA<-
  function(genExp, pheno, filename = "LIMMAresults"){
    "
    need pheno to be TP1, TP2, TP3, TP4, and ambi
    "
    
    require(limma)
    
    #subset genes to only informative ones
    sds<-apply(genExp,2,sd)
    if(length(which(sds==0))!=0){
      genExp<-genExp[,-which(sds==0)]
    }
    
    all_genes<-colnames(genExp)
    
    #LIMMA
    ### DEA setup: groups
    eset1 <- as.matrix(genExp[which(pheno!="Ambi"),])
    pheno[which(pheno!="Ambi")]-> Group
    dim(eset1)
    Group<-droplevels(Group)
    # eset2 <- (2^eset1)-1
    
    
    ###DEA
    design <- model.matrix(~0+Group)
    fit <- lmFit(t(eset1), design)
    contrast.matrix <-
      makeContrasts(
        "TP1" = GroupTP1 - (GroupTP2 + GroupTP3 + GroupTP4)/3,
        "TP2" = GroupTP2 - (GroupTP1 + GroupTP3 + GroupTP4)/3,
        "TP3" = GroupTP3 - (GroupTP1 + GroupTP2 + GroupTP4)/3,
        "TP4" = GroupTP4 - (GroupTP1 + GroupTP2 + GroupTP3)/3,
        levels = design
      )
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    ngenes = dim(eset1)[2]
    res = topTable(fit2,  n = ngenes)
    
    #write
    write.csv(res,paste0(filename,".csv"))
   
    res 
  }

fourgrpfgsea<-
  function(res,pathwaysname = "c7.all.v7.0", filename = paste0("fgsea_",pathwaysname)){
    #fgsea
    #http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2
    #MSigDB "C2: curated gene sets" 4762 gene sets
    #h.all.v7.0
    #c2.cp.reactome.v7.0
    #c2.cp.kegg.v7.0
    #c3.tft.v7.0
    #c3.mir.v7.0
    #c5.bp.v7.0
    #c5.cc.v7.0
    #c5.mf.v7.0
    #c7.all.v7.0
    
    require(data.table)
    require(fgsea)
    
    comps<-colnames(res)[1:4]
    
    #load gmt 
    gmt.file <- system.file("extdata", paste0(pathwaysname,".symbols.gmt"), package="fgsea")
    pathways <- gmtPathways(gmt.file)
    

    dir.create("fgsea_results")
    allres1<-c()
    #loop
    for(i in 1:length(comps)){
      # i=1
      my_ranks<-res[,i]
      names(my_ranks)<-rownames(res)
      
      fgseaRes <- fgsea(pathways = pathways, 
                        stats = my_ranks,
                        minSize=15,
                        maxSize=500,
                        nperm=10000)
      fgseaRes1 <- collapsePathways(fgseaRes,pathways = pathways, 
                                    stats = my_ranks)
      
      #write data table
      fwrite(fgseaRes, file=paste0("./fgsea_results/",comps[i],"_fgseaRes_",pathwaysname,".txt"), sep="\t", sep2=c("", " ", ""))
      write.csv(fgseaRes1$mainPathways, file=paste0("./fgsea_results/",comps[i],"_fgseaRes_collapse_",pathwaysname,".csv"))
      
      #make combined table
      colnames(fgseaRes)<-paste0(comps[i],"_",colnames(fgseaRes))
      allres1<-cbind.data.frame(allres1,fgseaRes)
    }
    
    
    # NES, padj, leadingEdge
    colnames(allres1)
    allres1<-data.frame(allres1)
    combres_ind<-c()
    for (i in 0:3){
      combres_ind<-c(combres_ind,c(i*8+c(3,1,6)+2))
    }
    subsetcombres<-allres1[,c(1,7,combres_ind)]
    colnames(subsetcombres)
    
    fwrite(allres1, file=paste0("./fgsea_results/",filename,"_allresults.txt"), 
           sep="\t", sep2=c("", " ", ""))
    fwrite(subsetcombres, file=paste0("./fgsea_results/",filename,"_shortresults.txt"), 
           sep="\t", sep2=c("", " ", ""))
    subsetcombres
  }

threegrpDEA<-
  function(genExp, pheno, filename = "LIMMAresults"){
    "
    need pheno to be tp1, tp2, neuroendocrine
    "
    
    require(limma)
    
    #subset genes to only informative ones
    sds<-apply(genExp,2,sd)
    if(length(which(sds==0))!=0){
      genExp<-genExp[,-which(sds==0)]
    }
    
    all_genes<-colnames(genExp)
    
    #LIMMA
    ### DEA setup: groups
    eset1 <- as.matrix(genExp)
    pheno-> Group
    # eset2 <- (2^eset1)-1
    
    
    ###DEA
    design <- model.matrix(~0+Group)
    fit <- lmFit(t(eset1), design)
    contrast.matrix <-
      makeContrasts(
        "tp1" = Grouptp1 - (Grouptp2 + Groupneuroendocrine)/2,
        "tp2" = Grouptp2 - (Grouptp1 + Groupneuroendocrine)/2,
        "neuroendocrine" = Groupneuroendocrine - (Grouptp1 + Grouptp2)/2,
        levels = design
      )
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    ngenes = dim(eset1)[2]
    res = topTable(fit2,  n = ngenes)
    
    #write
    write.csv(res,paste0(filename,".csv"))
    
    res 
  }


makeheatmap <- 
  function (cat1, cat2, db_umap_res, save = T, filename = "DBUheatmap"){
  my_sample_col<-data.frame(newpheno[,c("study2","DBUgrp")])
  annotation_colors<- list(
    LUAD_Group = c(`TP1`="#E41A1C",
                   `TP2`="#377EB8",
                   `TP3`="#4DAF4A",
                   `TP4`="#984EA3",
                   `Ambi`="Grey"),
    Study = c(`Dir`="#8DD3C7",
              `Lim`="#FFFFB3",
              `TCGA`="#BEBADA")
  )
  rownames(my_sample_col)<-rownames(newpheno)
  colnames(my_sample_col)[2]<-"LUAD_Group"
  colnames(my_sample_col)[1]<-"Study"
  
  #order
  my_sample_col1<-my_sample_col[order(my_sample_col$LUAD_Group),]
  dbumapgrps1<-db_umap_res[match(rownames(my_sample_col1),rownames(my_sample_col)),]
  rownames(dbumapgrps1)<-rownames(my_sample_col1)
  gaps<-as.numeric(table(my_sample_col1$LUAD_Group))
  for(i in 2:length(gaps)){gaps[i]<-gaps[i]+gaps[i-1]}
  
  if (save == T){
  png(paste0(filename,".png",width=6, height=10,units = "in",res=300))
  pheatmap(as.matrix(t(dbumapgrps1)),
           annotation_col = my_sample_col,
           annotation_colors = annotation_colors,
           cluster_rows = F,
           cluster_cols = F,
           gaps_col = gaps[1:4],
           annotation_names_row = F,
           annotation_names_col = F,
           fontsize_row = 1,
           fontsize_col = 0.5,
           color = RColorBrewer::brewer.pal(7, "BuPu")[2:7],
           cutree_cols = 4,
           show_rownames = F,
           show_colnames = F,
           clustering_distance_cols = "manhattan",
           width = 6,
           height = 5)
  dev.off()
  }
  else
    pheatmap(as.matrix(t(dbumapgrps1)),
             annotation_col = my_sample_col,
             annotation_colors = annotation_colors,
             cluster_rows = F,
             cluster_cols = F,
             gaps_col = gaps[1:4],
             annotation_names_row = F,
             annotation_names_col = F,
             fontsize_row = 1,
             fontsize_col = 0.5,
             color = RColorBrewer::brewer.pal(7, "BuPu")[2:6],
             cutree_cols = 4,
             show_rownames = F,
             show_colnames = F,
             clustering_distance_cols = "manhattan",
             width = 6,
             height = 5)
    
  }

plot_DBUiter<-
  function(umap_position_matrix, i, png = F, pdf = F){
    "gets the job done"

    x<-(i*2)-1
    df<-umap_position_matrix[,x:(x+1)]
    plot(df,xlab="",ylab="",pch=19,cex=0.7)
    
    if (png == T){
      png(paste0("02_TCGALUAD_DBUiteration#",i,"_unsupgrp_plot.png"),width = 3,height = 2.5,res = 300,units = "in")
      par(mar=c(2,2,0,0))
      plot(df, col=col1[samp_class],xlab="",ylab="",pch=19,cex=0.7)
      dev.off()
    }
    if (pdf == T){
      pdf(paste0("02_DBUiteration#",i,"_unsupgrp_plot.pdf"),width = 3,height = 2.5)
      par(mar=c(2,2,0,0))
      plot(df, col=col1[samp_class],xlab="",ylab="",pch=19,cex=0.7)
      dev.off()
    }
  }

plot_DBUiter2<-
  function(color, umap_position_matrix, i, png = T, pdf = F){
    "more fancy, needs more customization"
    color<-dplyr::recode_factor(color,
                         `TP1`=1,
                         `TP2`=2,
                         `TP3`=3,
                         `TP4`=4,
                         `Ambi`=5)
    samp_class<-as.numeric(color)
    col1<-RColorBrewer::brewer.pal(5,"Set1")
    col1[5]<-"Grey"
    
    x<-(i*2)-1
    df<-umap_position_matrix[,x:(x+1)]
    plot(df, col=col1[samp_class],xlab="",ylab="",pch=19,cex=0.7)
    
    if (png == T){
    png(paste0("02_TCGALUAD_DBUiteration#",i,"_unsupgrp_plot.png"),width = 3,height = 2.5,res = 300,units = "in")
    par(mar=c(2,2,0,0))
    plot(df, col=col1[samp_class],xlab="",ylab="",pch=19,cex=0.7)
    dev.off()
    }
    if (pdf == T){
    pdf(paste0("02_DBUiteration#",i,"_unsupgrp_plot.pdf"),width = 3,height = 2.5)
    par(mar=c(2,2,0,0))
    plot(df, col=col1[samp_class],xlab="",ylab="",pch=19,cex=0.7)
    dev.off()
    }
  }

#read all excel sheets function
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
#confidence function
conf_score <- function (a_num_vec){
  a_num_vec <- as.numeric(a_num_vec)
  a_num_vec <- sort(a_num_vec,decreasing = T)
  conf_score <- a_num_vec[1]/a_num_vec[2]
  # if(conf_score==Inf){conf_score = 100}
  return(conf_score)
}
#caller function
caller <- function (conf_score, named_freq_vec, threshold){
  if (conf_score>=threshold|conf_score==Inf){
    named_freq_vec <-sort(named_freq_vec,decreasing = T)
    caller = names(named_freq_vec)[1]
  } else if (conf_score<threshold){  
    named_freq_vec <-sort(named_freq_vec,decreasing = T)
    caller = paste0(names(named_freq_vec)[1],".",names(named_freq_vec)[2])
  }
  return(caller)
}

