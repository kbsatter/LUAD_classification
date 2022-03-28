


rm(list = ls())

# lung data tcga compare
library(tidyverse)
library(data.table)
library(magrittr)
library(umap)
library(factoextra)
library(pheatmap)
library(dbscan)



## getwd
df = fread("TCGA.LUNG.sampleMap_HiSeqV2.gz")
pheno = fread("TCGA.LUNG.sampleMap_LUNG_clinicalMatrix")
pheno %<>% filter(sample_type != "Solid Tissue Normal")
pheno %<>% mutate(Tumor.type = ifelse(str_detect(histological_type, "Adeno"), 
                               "Adenocarcinoma", "Other"))


df2 = df[ , colnames(df) %chin% pheno$sampleID, with = F]
pheno2 = pheno[ match(colnames(df2), pheno$sampleID), ]
identical(pheno2$sampleID, colnames(df2))

GeneID = df$sample
Histology = pheno2$Tumor.type
u1 = umap(t(df2))
u1$layout %>% as.data.frame %>% 
  ggplot(aes(V1,V2, color = Histology)) +
  geom_point(alpha = 0.5) + theme_bw() +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(text =element_text(size = 15))




####
rm(list = ls())

# lung data tcga compare
library(tidyverse)
library(magrittr)
library(umap)
library(factoextra)
library(pheatmap)
library(dbscan)



## getwd
df = read.csv("TCGA.LUNG.sampleMap_HiSeqV2.gz", sep = "\t", row.names = 1)
pheno = read.csv("TCGA.LUNG.sampleMap_LUNG_clinicalMatrix", sep = '\t')
pheno %<>% filter(sample_type != "Solid Tissue Normal")
pheno$sampleID = gsub("-", ".", pheno$sampleID)
df = df[ , colnames(df) %in% pheno$sampleID]
pheno = pheno[ match(colnames(df), pheno$sampleID), ]
sdata = read.csv("survival_LUNG_survival.txt", sep = '\t')
sdata = sdata[ match(pheno2$X_PATIENT, sdata$X_PATIENT), ]

### 
### unsupervised
source("./Lung/scripts/myfunctions.R")
colorlabel = pheno$histological_type
df = t(df)
for (i in c("euclidean", "manhattan", "cosine")){
  #i="euclidean"
  pdf(paste0("01b_Umap_parameters_grid_search_",i,"_",gsub("-","",Sys.Date()),".pdf"))
  umap_parameters_gridsearch(df,colnames(df),colorlabel)
  dev.off()
}


###
umap_position_matrix2 = umapiterations(dataset = df,
                                       genes = colnames(df),
                                       no_iters = 1000,
                                       rand_gene_no = 1000,
                                       n_neighbors = 50,
                                       min_dist = 0.1,
                                       metric = "manhattan")



Apar(mfrow = c(3,3))
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
k = 2
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

pheno = cbind(pheno, consens_calls)
table(pheno$histological_type, pheno$samp_clust6)

pheno %<>% mutate(Histology = ifelse(str_detect(histological_type, "Adeno"), 
                               "Adenocarcinoma", "Other"))
table(pheno$Histology, pheno$samp_clust6)
table(pheno$sample_type)

write.csv(pheno, file = "Phenotype_all_lungcancer.csv")

save(db_umap_res, pheno, file = "alllung2.Rdata")

library(ggalluvial)
pheno %>% 
  ggplot( aes(axis1 = Histology, axis3 = final_cluster)) +
  scale_x_discrete(limits = c("Histology", "DBU cluster"))+
  xlab("") +
  geom_alluvium(aes(fill = Histology)) +
  geom_stratum() + theme_classic() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + theme(text = element_text(size = 15, color = "black")) 
