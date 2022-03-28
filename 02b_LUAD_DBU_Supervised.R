

# Analysis Date: June 9, 2020
# Author: Paul Tran
# Title: Density based Clustering of UMAP iterations for Lung Cancer

"
TO DO: 
get glmnet variables and gene names
plot alluvial

load combined data

combined three datasets for group discovery more powered to find subgroups.
low likelihood of one driving since all three are similar size
if only one dataset has samples in that group, more likely artifact, otherwise may be real

run UMAP
run DBSCAN
perform consensus clustering and make group calls
plot representative umap and heatmap
save
"

rm(list=ls(all=TRUE))

#set seed
set.seed(123)

#set directory
my.dir<-"/Users/paultran/Box/TCGA Cancer Classification Project/Lung 20200420"
setwd(my.dir)

#load functions
source("./scripts/myfunctions.R")

#load libraries
library(umap)
library(edgeR)
library(dplyr)
library(caret)
library(doParallel)
require(openxlsx)

##############  1) load and process data #######
load("./data/Merged_genexp_pheno.Rdata")
TCGA_exp<-genExp_scaled[pheno_merged$Study=="TCGA",]#509
Dir_exp<-genExp_scaled[pheno_merged$Study=="Director",]#443
Lim_exp<-genExp_scaled[pheno_merged$Study=="Lim",]#632


############### 2) TCGA UMAP ##################

#2a. identify optimal umap parameters for each study
for(i in c("TCGA","Director","Lim")){
pdf(paste0("./Results/",i,"_umap_gridsearch.pdf"))
umap_parameters_gridsearch(dataset = genExp_scaled[pheno_merged$Study==i,],
                           genelist = colnames(genExp_scaled),
                           colorlabel = 1)
dev.off()
}

#2b. Run 1000 umap iterations for TCGA
umap_position_matrix2<-umapiterations(dataset = TCGA_exp, 
                                      genes = colnames(TCGA_exp),
                                      no_iters = 1000,
                                      rand_gene_no = 5000,
                                      n_neighbors = 20,
                                      min_dist = 0.1,
                                      metric = "manhattan")

#plot to test
par(mfrow = c(3,3))
for(i in 1:9){
  plot_DBUiter(umap_position_matrix2,i)
}


############### 3) TCGA DBSCAN ##################
#3c. optimize DBSCAN parameters

dbutest(umap_position_matrix = umap_position_matrix2,
        iter = 70,
        k = 20,
        eps = 0.65)
#3d. run DBSCAN on all umap iters
db_umap_res<-dbuiters(umap_position_matrix = umap_position_matrix2,
                      k = 20,
                      eps = 0.65)

############### 4) TCGA Consensus clustering ##################
#Find pairwise co-grouping frequencies for each sample
Mfreq<-group_consensus_freq(umap_position_matrix2,db_umap_res)
#make group calls based on clustering pairwise frequencies
k=5
pheatmap::pheatmap(Mfreq,scale = "none",
                   clustering_distance_cols = "manhattan",
                   clustering_distance_rows = "manhattan",
                   cutree_rows = k,
                   cutree_cols = k)

consens_calls<-makegrpcalls_clust(Mfreq,k)
table(consens_calls$samp_clust6)

#renumber groups
consens_calls$final_cluster<-change_grp_names4(consens_calls$samp_clust6)
consens_calls$final_cluster[consens_calls$samp_clus_div_res>7]<-NA
table(consens_calls$final_cluster)

#merge with pheno
x=5
plot(umap_position_matrix2[,x:(x+1)], col=consens_calls$final_cluster,xlab="",ylab="",pch=19,cex=0.7)

############### 5) TCGA Supervised Classication of DBU ##################
#train data and holdout
TCGA_exp %>% filter(!is.na(consens_calls$final_cluster))->TCGA_exp_grouped
TCGA_exp_grouped$class<-recode_factor(consens_calls$final_cluster[!is.na(consens_calls$final_cluster)],
              `1` = "TP1",
              `2` = "TP2",
              `3` = "TP3",
              `4` = "TP4")

TCGA_idx = createDataPartition(TCGA_exp_grouped$class, p = 0.75, list = FALSE)
TCGA_trn = TCGA_exp_grouped[TCGA_idx, ] #377
TCGA_tst = TCGA_exp_grouped[-TCGA_idx, ] #123

#glmnet
set.seed(123)
# set up parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)
getDoParWorkers()

##########pam##########
pam_mod1 = train(
  class~.,
  data = TCGA_trn,
  method = "pam",
  tuneGrid = expand.grid(threshold = c(0.01,0.1,1,5,10)),
  trControl = trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = multiClassSummary)
)

#holdout validation
mean(predict(pam_mod1,TCGA_tst)==TCGA_tst$class) #Accuracy 0.89
confusionMatrix(data=predict(pam_mod1,TCGA_tst),reference = as.factor(TCGA_tst$class))

#glmnet model
my_glmnet_model <- train(
  class ~ .,
  method = "glmnet",
  data = TCGA_trn,
  trControl = trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = multiClassSummary)
)

#getting the coefficients of the final model
coefficients <- coef(my_glmnet_model$finalModel, my_glmnet_model$bestTune$lambda)
str(coefficients$TP1) #390
str(coefficients$TP2) #263
str(coefficients$TP3) #248
str(coefficients$TP4) #194


#holdout validation
mean(predict(my_glmnet_model,TCGA_tst)==TCGA_tst$class) #Accuracy 0.94
confusionMatrix(data=predict(my_glmnet_model,TCGA_tst),reference = as.factor(TCGA_tst$class))

#get accuray results
results2<-resamples(list(pam_mod=pam_mod1,glmnet_mod=my_glmnet_model))
summary(results2);bwplot(results2)
mean(results2$values$`glmnet_mod~Accuracy`)
norm.interval(results2$values$`glmnet_mod~Accuracy`)
mean(results2$values$`pam_mod~Accuracy`)
norm.interval(results2$values$`pam_mod~Accuracy`)

#plot accuracies
ggplot(data = data.frame(Accuracy=results2$values$`glmnet_mod~Accuracy`,
                         Glmnet = rep("Glmnet",50)),
       aes(y=Accuracy,x=Glmnet))+
  geom_boxplot()+
  geom_jitter()+
  geom_point(aes(x="Glmnet",y=0.99),colour="red")+
  coord_cartesian(ylim = c(0,1))+
  theme_classic()
ggsave("./Results/03_DBUsup_acc_bxplt.png",width = 2.5,height = 2.5)

#test set 1
multi_ROC_4grp(as.vector(TCGA_tst$class), my_glmnet_model,TCGA_tst[,-6512]) # 0.987, 0.989, 0.999, 0.985, 0.988, 0.987
plot_multi_ROC_4grp(as.vector(TCGA_tst$class), my_glmnet_model,TCGA_tst[,-6512])
ggsave("./Results/02_TCGAtst_ROC.png",width = 2.5,height = 2.5)

multi_PR_4grp(as.vector(TCGA_tst$class), my_glmnet_model,TCGA_tst[,-6512]) #0.996,0.811, 0.992, 0.893
plot_multi_PR_4grp(as.vector(TCGA_tst$class), my_glmnet_model,TCGA_tst[,-6512])
ggsave("./Results/02_TCGAtst_pr.png",width = 2.5,height = 2.5)

#classify all
pheno_merged$TCGA_sup_class<-predict(my_glmnet_model,genExp_scaled)
table(pheno_merged$TCGA_sup_class)


############### 6) TCGA SAVE ##################
write.csv(umap_position_matrix2[1:dim(umap_position_matrix2)[1],
                                1:dim(umap_position_matrix2)[2]],
          "./data/TCGA_umap_posmat.csv")
save(db_umap_res, consens_calls, file = "./data/TCGA_dbu_grpmat.Rdata")
saveRDS(my_glmnet_model,"./data/TCGA_glmnet_model.rds")

list_of_datasets <- list("TCGA1_glmnet_coef" = cbind(colnames(TCGA_exp)[coefficients$TP1@i+1],coefficients$TP1@x), 
                         "TCGA2_glmnet_coef" = cbind(colnames(TCGA_exp)[coefficients$TP2@i+1],coefficients$TP2@x),
                         "TCGA3_glmnet_coef" = cbind(colnames(TCGA_exp)[coefficients$TP3@i+1],coefficients$TP3@x),
                         "TCGA4_glmnet_coef" = cbind(colnames(TCGA_exp)[coefficients$TP4@i+1],coefficients$TP4@x))
write.xlsx(list_of_datasets, file = "./results/TCGA_glmnet_coefficients.xlsx")

################### Dir ######################
################### 7) Dir UMAP ##############
#2b. Run 1000 umap iterations for TCGA
umap_position_matrix2<-umapiterations(dataset = Dir_exp, 
                                      genes = colnames(Dir_exp),
                                      no_iters = 1000,
                                      rand_gene_no = 5000,
                                      n_neighbors = 20,
                                      min_dist = 0.1,
                                      metric = "manhattan")

#plot to test
par(mfrow = c(3,3))
for(i in 1:9){
  plot_DBUiter(umap_position_matrix2,i)
}


############### 3) Dir DBSCAN ##################
#3c. optimize DBSCAN parameters

dbutest(umap_position_matrix = umap_position_matrix2,
        iter = 700,
        k = 20,
        eps = 0.6)
#3d. run DBSCAN on all umap iters
db_umap_res<-dbuiters(umap_position_matrix = umap_position_matrix2,
                      k = 20,
                      eps = 0.6)

############### 4) Dir Consensus clustering ##################
#Find pairwise co-grouping frequencies for each sample
Mfreq<-group_consensus_freq(umap_position_matrix2,db_umap_res)
#make group calls based on clustering pairwise frequencies
k=6
pheatmap::pheatmap(Mfreq,scale = "none",
                   clustering_distance_cols = "manhattan",
                   clustering_distance_rows = "manhattan",
                   cutree_rows = k,
                   cutree_cols = k)

consens_calls<-makegrpcalls_clust(Mfreq,k)
table(consens_calls$samp_clust6)

#renumber groups
consens_calls$final_cluster<-change_grp_names4(consens_calls$samp_clust6)
consens_calls$final_cluster[consens_calls$samp_clus_div_res>10]<-NA
table(consens_calls$final_cluster)

#merge with pheno
x=5
plot(umap_position_matrix2[,x:(x+1)], col=consens_calls$final_cluster,xlab="",ylab="",pch=19,cex=0.7)

############### 5) Dir Supervised Classication of DBU ##################
#train data and holdout
Dir_exp %>% filter(!is.na(consens_calls$final_cluster))->Dir_exp_grouped
Dir_exp_grouped$class<-recode_factor(consens_calls$final_cluster[!is.na(consens_calls$final_cluster)],
                                      `1` = "TP1",
                                      `2` = "TP2",
                                      `3` = "TP3",
                                      `4` = "TP4")
set.seed(123)
Dir_idx = createDataPartition(Dir_exp_grouped$class, p = 0.75, list = FALSE)
Dir_trn = Dir_exp_grouped[Dir_idx, ] #369
Dir_tst = Dir_exp_grouped[-Dir_idx, ] #122

#glmnet
# set up parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)
getDoParWorkers()

##########pam##########
pam_mod1 = train(
  class~.,
  data = Dir_trn,
  method = "pam",
  tuneGrid = expand.grid(threshold = c(0.01,0.1,1,5,10)),
  trControl = trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = multiClassSummary)
)

#holdout validation
mean(predict(pam_mod1,Dir_tst)==Dir_tst$class) #Accuracy 0.95
confusionMatrix(data=predict(pam_mod1,Dir_tst),reference = as.factor(Dir_tst$class))

#glmnet model
my_glmnet_model <- train(
  class ~ .,
  method = "glmnet",
  data = Dir_trn,
  trControl = trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 5,
                           summaryFunction = multiClassSummary)
)

#getting the coefficients of the final model
coefficients <- coef(my_glmnet_model$finalModel, my_glmnet_model$bestTune$lambda)
str(coefficients$TP1) #486
str(coefficients$TP2) #315
str(coefficients$TP3) #274
str(coefficients$TP4) #337

#holdout validation
mean(predict(my_glmnet_model,Dir_tst)==Dir_tst$class) #Accuracy 0.96
confusionMatrix(data=predict(my_glmnet_model,Dir_tst),reference = as.factor(Dir_tst$class))

#get accuray results
results2<-resamples(list(pam_mod=pam_mod1,glmnet_mod=my_glmnet_model))
summary(results2);bwplot(results2)
mean(results2$values$`glmnet_mod~Accuracy`)
norm.interval(results2$values$`glmnet_mod~Accuracy`)
mean(results2$values$`pam_mod~Accuracy`)
norm.interval(results2$values$`pam_mod~Accuracy`)

#plot accuracies
ggplot(data = data.frame(Accuracy=results2$values$`glmnet_mod~Accuracy`,
                         Glmnet = rep("Glmnet",50)),
       aes(y=Accuracy,x=Glmnet))+
  geom_boxplot()+
  geom_jitter()+
  geom_point(aes(x="Glmnet",y=0.99),colour="red")+
  coord_cartesian(ylim = c(0,1))+
  theme_classic()
ggsave("./Results/03_DBUsup_acc_bxplt.png",width = 2.5,height = 2.5)

#test set 1
multi_ROC_4grp(as.vector(Dir_tst$class), my_glmnet_model,Dir_tst[,-6512]) # 0.997, 0.997, 0.9999, 0.998, 0.997, 0.997

plot_multi_ROC_4grp(as.vector(Dir_tst$class), my_glmnet_model,Dir_tst[,-6512])
ggsave("./Results/02_Dirtst_roc.png",width = 2.5,height = 2.5)

multi_PR_4grp(as.vector(Dir_tst$class), my_glmnet_model,Dir_tst[,-6512]) #0.9996, 0.86, 0.998, 0.90, 0.93, 0.99
plot_multi_PR_4grp(as.vector(Dir_tst$class), my_glmnet_model,Dir_tst[,-6512])
ggsave("./Results/02_Dirtst_pr.png",width = 2.5,height = 2.5)

#classify all
pheno_merged$Dir_sup_class<-predict(my_glmnet_model,genExp_scaled)
table(pheno_merged$Dir_sup_class, pheno_merged$TCGA_sup_class)


############### 5) Dir SAVE ##################
write.csv(umap_position_matrix2[1:dim(umap_position_matrix2)[1],
                                1:dim(umap_position_matrix2)[2]],
          "./data/Dir_umap_posmat.csv")
save(db_umap_res, consens_calls, file = "./data/Dir_dbu_grpmat.Rdata")
saveRDS(my_glmnet_model,"./data/Dir_glmnet_model.rds")

list_of_datasets <- list("Dir1_glmnet_coef" = cbind(colnames(Dir_exp)[coefficients$TP1@i+1],coefficients$TP1@x), 
                         "Dir2_glmnet_coef" = cbind(colnames(Dir_exp)[coefficients$TP2@i+1],coefficients$TP2@x),
                         "Dir3_glmnet_coef" = cbind(colnames(Dir_exp)[coefficients$TP3@i+1],coefficients$TP3@x),
                         "Dir4_glmnet_coef" = cbind(colnames(Dir_exp)[coefficients$TP4@i+1],coefficients$TP4@x))
write.xlsx(list_of_datasets, file = "./results/Dir_glmnet_coefficients.xlsx")


############### 
############# 2) classify samples using Wilkerson et al ####
#http://cancer.unc.edu/nhayes/publications/adenocarcinoma.2012/
load("./data/TCGA_Lim_Dir_genexp_pheno.Rdata")
wilk_centroids<-read.csv("./data/wilkerson.2012.LAD.predictor.centroids.csv",row.names = 1)
common_genes<-intersect(colnames(genExp_scaled),rownames(wilk_centroids))
TCGA_wilk_exp<-genExp_scaled[pheno_merged$Study=="TCGA"&!is.na(TCGA_pheno$paper_expression_subtype),]
table(TCGA_pheno$paper_expression_subtype)

wilk_centroids2<-wilk_centroids[common_genes,]#447

wilk_cor_res<-cbind.data.frame(grp=TCGA_pheno$paper_expression_subtype[!is.na(TCGA_pheno$paper_expression_subtype)],
                               tru=apply(TCGA_wilk_exp[,common_genes],1,function(x)cor(x,wilk_centroids2$bronchioid)),
                               pp=apply(TCGA_wilk_exp[,common_genes],1,function(x)cor(x,wilk_centroids2$magnoid)),
                               pi=apply(TCGA_wilk_exp[,common_genes],1,function(x)cor(x,wilk_centroids2$squamoid)))

set.seed(430)
wilk_idx = createDataPartition(wilk_cor_res$grp, p = 0.75, list = FALSE)
wilk_trn = wilk_cor_res[wilk_idx, ] #172
wilk_tst = wilk_cor_res[-wilk_idx, ] #55

multinom_wilk_mod = train(
  grp~.,
  data = wilk_trn,
  method = "multinom",
  trControl = trainControl(method = "repeatedcv",number = 10, repeats = 5,summaryFunction = multiClassSummary)
)

glmnet_wilk_mod = train(
  grp~.,
  data = wilk_trn,
  method = "glmnet",
  trControl = trainControl(method = "repeatedcv",number = 10, repeats = 5,summaryFunction = multiClassSummary)
)

#get accuray results
results<-resamples(list(multinom=multinom_wilk_mod,glmnet=glmnet_wilk_mod))
summary(results)
mean(results$values$`multinom~Accuracy`)
norm.interval(results$values$`multinom~Accuracy`)
mean(results$values$`glmnet~Accuracy`)
norm.interval(results$values$`glmnet~Accuracy`)

mean(predict(multinom_wilk_mod,wilk_tst[,-1])==wilk_tst$grp) #Accuracy 0.96
confusionMatrix(data=predict(multinom_wilk_mod,wilk_tst[,-1]),reference = wilk_tst$grp)

mean(predict(glmnet_wilk_mod,wilk_tst[,-1])==wilk_tst$grp) #Accuracy 0.96
confusionMatrix(data=predict(glmnet_wilk_mod,wilk_tst[,-1]),reference = wilk_tst$grp)

#plot accuracies
ggplot(data = data.frame(Accuracy=results$values$`multinom~Accuracy`,
                         Multinom = rep("Multinom",50)),
       aes(y=Accuracy,x=Multinom))+
  geom_boxplot()+
  geom_jitter()+
  geom_point(aes(x="Multinom",y=0.96),colour="red")+
  coord_cartesian(ylim = c(0,1))+
  theme_classic()
ggsave("./Results/03_Wilkerson_acc_bxplt.png",width = 2.5,height = 2.5)

#plot val roc
multi_ROC_3grp(as.vector(wilk_tst$grp), multinom_wilk_mod,wilk_tst[,-1]) # 0.997, 0.997, 0.9999, 0.998, 0.997, 0.997

plot_multi_ROC_3grp(as.vector(wilk_tst$class), multinom_wilk_mod,wilk_tst[,-4])
ggsave("./Results/02_Dirtst_roc.png",width = 2.5,height = 2.5)

multi_PR_3grp(as.vector(wilk_tst$class), multinom_wilk_mod,wilk_tst[,-6512]) #0.9996, 0.86, 0.998, 0.90, 0.93, 0.99
plot_multi_PR_3grp(as.vector(wilk_tst$class), multinom_wilk_mod,wilk_tst[,-6512])
ggsave("./Results/02_Dirtst_pr.png",width = 2.5,height = 2.5)


roc_data<-data.frame(
  G1_true=fastDummies::dummy_cols(wilk_tst$grp)[,2],
  G2_true=fastDummies::dummy_cols(wilk_tst$grp)[,3],
  G3_true=fastDummies::dummy_cols(wilk_tst$grp)[,4],
  G1_pred_m1=predict(multinom_wilk_mod,wilk_tst[,-1],type = "prob")$`prox.-inflam`,
  G2_pred_m1=predict(multinom_wilk_mod,wilk_tst[,-1],type = "prob")$`prox.-prolif.`,
  G3_pred_m1=predict(multinom_wilk_mod,wilk_tst[,-1],type = "prob")$TRU)

roc_test <- multiROC::multi_roc(roc_data) #0.99, 0.99, 0.997, 0.99, 0.99
plot_roc_df<-multiROC::plot_roc_data(roc_test)
torm<-which(plot_roc_df$Group=="Micro"|plot_roc_df$Group=="Macro")
plot_roc_df$Group<-dplyr::recode(plot_roc_df$Group, G1 = "prox.-inflam", G2 = "prox.-prolif.", G3 = "TRU")

g<-ggplot(plot_roc_df[-torm,], aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
g
ggsave("./Results/02_Wilkerson_ROC.png",
       width = 2.5,
       height = 2.5,
       dpi = 300,
       units = "in")

#run on all
all_cor_res<-data.frame(tru = apply(genExp_scaled[,common_genes],1,function(x)cor(x,wilk_centroids2$bronchioid)),
                        pp = apply(genExp_scaled[,common_genes],1,function(x)cor(x,wilk_centroids2$magnoid)),
                        pi = apply(genExp_scaled[,common_genes],1,function(x)cor(x,wilk_centroids2$squamoid)))

all_cor_res$wilk_class<-predict(multinom_mod,all_cor_res)

pheno_merged<-cbind.data.frame(pheno_merged,all_cor_res)


################# alluvial plot #########3
clinfactors_forsankey<-data.frame(
  table(
    pheno_merged$DBU_sup_class,
    pheno_merged$wilk_class
  )
)

colnames(clinfactors_forsankey) <-
  c("DBU_Class",
    "Wilkerson_Class",
    "Freq"
  )
clinfactors_forsankey$Wilkerson_Class = 
  factor(clinfactors_forsankey$Wilkerson_Class,
         levels(clinfactors_forsankey$Wilkerson_Class)[c(1,3,2)])


png("./Results/03_DBU_wilkerson_alluvial.png", width = 6, height = 3, res = 300,units = "in")
alluvial::alluvial(clinfactors_forsankey[,1:2], freq=clinfactors_forsankey$Freq,
                   hide = clinfactors_forsankey$Freq < 5,
                   col = c("#E41A1C",
                           "#377EB8",
                           "#4DAF4A",
                           "#984EA3"),
                   blocks=T,cw=0.25,cex=0.5
)
dev.off()




#### 4) representative umap plot and heatmap #############
#umap
png(paste0("./Results/02_LUAD_umapplot.png"),width = 3,height = 2.5,res = 300,units = "in")
par(mar=c(2,2,0,0))

colors<-recode(consens_calls$final_cluster,
       `1` = "#E41A1C",
       `2` = "#377EB8",
       `3` = "#4DAF4A",
       `4` = "#984EA3")
colors[is.na(colors)]<-"gray90"

plot(umap_position_matrix2[,1:2],col=scales::alpha(colors,0.8),pch=19,cex=0.4)
dev.off()



#plot heatmap
pheno_merged<-cbind.data.frame(pheno_merged,consens_calls)

my_samp_col<-data.frame(Old_Cluster = paste0("Group_",pheno_merged$samp_clust6),
                        Cluster_Dispersion = pheno_merged$samp_clus_div_res,
                        Cluster = paste0("LUAD_",ifelse(is.na(pheno_merged$final_cluster),"Ambi",pheno_merged$final_cluster)),
                        Study = pheno_merged$Study
                        )
rownames(my_samp_col)<-pheno_merged$samp_id
rownames(db_umap_res3_4)<-pheno_merged$samp_id

my_iter_col<-data.frame(Cluster_No = c(rep("3",dim(db_umap_res3)[2]),rep("4",dim(db_umap_res4)[2])))
rownames(my_iter_col)<-colnames(db_umap_res3_4)

my_colors<-list(
  Cluster = c(LUAD_1 = "#E41A1C",
              LUAD_2 = "#377EB8",
              LUAD_3 = "#4DAF4A",
              LUAD_4 = "#984EA3",
              LUAD_Ambi = "gray90")
)

h<-pheatmap::pheatmap(db_umap_res3_4, 
                      scale = "none",
                      clustering_distance_rows = "manhattan",
                      cutree_rows = 6,
                      annotation_col = my_iter_col,
                      annotation_row = my_samp_col,
                      annotation_colors = my_colors,
                      show_rownames = F,
                      show_colnames = F,
                      )

ggplot2::ggsave("./Results/TCGA_Lim_Dir_dbu_heatmap.png",plot = h)


#### save final pheno
save(pheno_merged, file = "./data/TCGA_Lim_Dir_pheno.Rdata")

