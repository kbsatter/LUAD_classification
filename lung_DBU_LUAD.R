


rm(list = ls())

# lung data tcga compare
library(tidyverse)
library(magrittr)
library(umap)
library(factoextra)
library(pheatmap)
library(dbscan)



## getwd
df = read.csv("TCGA.LUAD.sampleMap_HiSeqV2.gz", sep = "\t", row.names = 1)
pheno = read.csv("Phenotype_LUAD-VS-LUSC.csv")
pheno %<>% filter(sample_type != "Solid Tissue Normal")
pheno2 = pheno %>% filter(Histology == "Adenocarcinoma" &
                    final_cluster == 4)


df2 = df[ , colnames(df) %in% pheno2$sampleID]
pheno2 = pheno2[ match(colnames(df2), pheno2$sampleID), ]


sdata = read.csv("survival_LUNG_survival.txt", sep = '\t')
sdata = sdata[ match(pheno2$X_PATIENT, sdata$X_PATIENT), ]

### 
### unsupervised
source("./Lung/scripts/myfunctions.R")
colorlabel = pheno2$histological_type
df2 = t(df2)
for (i in c("euclidean", "manhattan", "cosine")){
  #i="euclidean"
  pdf(paste0("01b_Umap_parameters_grid_search_",i,"_",gsub("-","",Sys.Date()),".pdf"))
  umap_parameters_gridsearch(df,colnames(df),colorlabel)
  dev.off()
}


###
umap_position_matrix2 = umapiterations(dataset = df2,
                                       genes = colnames(df2),
                                       no_iters = 1000,
                                       rand_gene_no = 5000,
                                       n_neighbors = 20,
                                       min_dist = 0.1,
                                       metric = "cosine")



par(mfrow = c(3,3))
for(i in 1:9){
  plot_DBUiter(umap_position_matrix2,i)
}

dev.off()
#3c. optimize DBSCAN parameters

dbutest(umap_position_matrix = umap_position_matrix2,
        iter = 50,
        k = 5,
        eps = 0.4)

db_umap_res<-dbuiters(umap_position_matrix = umap_position_matrix2,
                      k = 5,
                      eps = 0.4)
Mfreq <- group_consensus_freq(umap_position_matrix2,db_umap_res)
#make group calls based on clustering pairwise frequencies
k = 5
pheatmap::pheatmap(Mfreq,scale = "none",
                   clustering_distance_cols = "manhattan",
                   clustering_distance_rows = "manhattan",
                   cutree_rows = k,
                   cutree_cols = k
                   )

consens_calls<-makegrpcalls_clust(Mfreq,k)
table(consens_calls$samp_clust6)

consens_calls$final_cluster<-change_grp_names5(consens_calls$samp_clust6)
consens_calls$final_cluster[consens_calls$samp_clus_div_res>7]<-NA
table(consens_calls$final_cluster)
# 
x = 50
plot(umap_position_matrix2[,x:(x+1)], col=consens_calls$final_cluster,
     xlab="",ylab="",pch=19,cex=0.7,
     main = "UMAP iteration 50")
plot(umap_position_matrix2[,400:(400+1)], col=consens_calls$final_cluster,
     xlab="",ylab="",pch=19,cex=0.7,
     main = "UMAP iteration 400")
plot(umap_position_matrix2[,750:(750+1)], col=consens_calls$final_cluster,
     xlab="",ylab="",pch=19,cex=0.7,
     main = "UMAP iteration 750")

colnames(consens_calls) = c("init_call", "cl_divergence", "LUAD_SG")

pheno2 = cbind(pheno2, consens_calls)
table(pheno2$Histology, pheno2$LUAD_SG)

write.csv(pheno2, "LUAD_Subgroup.csv")
save(db_umap_res, pheno2, file = "LUAD_SG_UMAP_res.Rdata")


### survival
library(survival)
library(survminer)
identical(sdata$X_PATIENT, pheno2$X_PATIENT)

sdata = cbind.data.frame(sdata, consens_calls)
sdata$LUAD_SG = as.factor(sdata$LUAD_SG)
head(sdata)
sdata %<>% mutate(OS.T = OS.time/30)
sdata %<>% filter(LUAD_SG != 5)
sdata$OS.T[sdata$OS.T > 60] = 60
sdata = sdata[!is.na(sdata$OS.T), ]
msurv = survfit(Surv(OS.T, OS) ~ LUAD_SG, data = sdata)
ggsurvplot(msurv, linetype = "strata", #legend = c(0.2,0.2), 
           size = 1.4, risk.table = T, legend.labs = c("TP1", "TP2", "TP3", "TP4"),
           pval = T)

sdata2 = merge(sdata, pheno2[ , c( "X_PATIENT", "pathologic_stage")], 
               all.x = T, all.y = F, by = "X_PATIENT")
msurv2 = survfit(Surv(OS.T, OS) ~ pathologic_stage, data = sdata2)

table(sdata2$pathologic_stage)
sdata2 %<>% mutate(Stage = ifelse(pathologic_stage == "Stage I" | 
                                   pathologic_stage == "Stage IA" |
                                   pathologic_stage == "Stage IB", "Stage I",
                  ifelse(pathologic_stage == "Stage II" |
                           pathologic_stage == "Stage IIA" |
                           pathologic_stage == "Stage IIB", "Stage II", 
                  ifelse(pathologic_stage == "Stage IIIA" |
                           pathologic_stage == "Stage IIIB", "Stage III", 
                    ifelse(pathologic_stage == "Stage IV", "Stage IV", "Other")))))
sdata2 %<>% filter(Stage != "Other")
msurv2 = survfit(Surv(OS.T, OS) ~ Stage, data = sdata2)
ggsurvplot(msurv2)
