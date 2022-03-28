
# Analysis Date: May 12, 2020
# Author: Paul Tran
# Title: Density based Clustering of UMAP iterations for Lung Cancer

"
fig 1
load, clean, filter, standardize TCGA
  can't standardize by gene bc it leave na
load,clean, filter, standardize lim
load, clean, filter, standardize dir

combined three datasets for group discovery more powered to find subgroups.
low likelihood of one driving since all three are similar size
if only one dataset has samples in that group, more likely artifact, otherwise may be real

combine pheno files

plot pca and dendrogram showing good merge

save
"

rm(list=ls(all=TRUE))

#set seed
set.seed(123)

#set directory
my.dir<-"/Users/paultran/Box/TCGA Cancer Classification Project/Lung 20200420"
setwd(my.dir)

#load functions
source("myfunctions.R")

#load libraries
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(umap)
library(GEOquery)
require(Biobase)
library(edgeR)
library(dplyr)


############## 1) Download/load and process data #######
#############  1A) TCGA         #######################
"
https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
https://rdrr.io/bioc/GDCRNATools/src/R/gdcVoomNormalization.R
"

query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Primary Tumor")

GDCdownload(query, 
            directory = paste0(my.dir,"/data"))

data <- GDCprepare(query, 
                   directory = paste0(my.dir,"/data"),
                   save = T,
                   save.filename = "./data/TCGA_LUAD.RData")

#remove ffpe and duplicated samples
TCGA_genExp<-assay(data)[,!colData(data)$is_ffpe&!duplicated(colData(data)$patient)]
rownames(TCGA_genExp)<-rowRanges(data)$external_gene_name

#make TCGA pheno file
TCGA_pheno<-colData(data)[!colData(data)$is_ffpe&!duplicated(colData(data)$patient),]


#tmm and voom
expr = DGEList(counts = TCGA_genExp)
expr = calcNormFactors(expr)
  
## filter out low expression genes
keepALL <- rowSums(cpm(expr) > 1) >= 0.5*ncol(TCGA_genExp)
exprALL <- expr[keepALL,,keep.lib.sizes = TRUE]

TCGA_genExp <- voom(exprALL, design=NULL, plot = FALSE)$E
  

############# 1B) LIM      ########################
Lim_genExp<-read.delim("./data/Lim_expr.txt",row.names = 1)
Lim_pheno<-read.delim("./data/Lim_pheno.txt",row.names = 1)

#resort genExp
Lim_genExp<-Lim_genExp[,rownames(Lim_pheno)]

#keep only ADC
which(Lim_pheno$Histology..ADC..adenocarcinoma..LCC..large.cell.carcinoma..SCC..squamous.cell.carcinoma.=="ADC")->sub1_ind
Lim_pheno[sub1_ind,]->Lim_pheno1
Lim_genExp[,sub1_ind]->Lim_genExp1


############# 1C) DIRECTOR   ###############
## download the NCI Director's Challenge Dataset
gse <- getGEO("GSE68465",GSEMatrix=F)
pheno<-pData(getGEO("GSE68465",GSEMatrix=T)[[1]])

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id=='GPL96'},GSMList(gse))

# get the probeset ordering
probesets <- Table(GPLList(gse)[[1]])$ID

# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]

# make a compliant ExpressionSet
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno1 <- as(pdata,"AnnotatedDataFrame")
eset2 <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno1)

# gene names
gene_symb<-Table(GPLList(gse)[[1]])$`Gene Symbol`
mydat<-cbind.data.frame(data.matrix,gene=gene_symb)
mydat.max <- aggregate(. ~ gene, data = mydat, max)
rownames(mydat.max)<-mydat.max$gene
mydat.max<-mydat.max[,-1]

# make dir
Dir_genExp<-mydat.max[,pheno$`disease_state:ch1`=="Lung Adenocarcinoma"]
Dir_pheno<-pheno[pheno$`disease_state:ch1`=="Lung Adenocarcinoma",]



############## 2) UMAP and COMBINE DATA ###################
dir.create("./Results")
pdf("./Results/UMAP_TCGA_Lim_Dir_comb_scaled.pdf")
par(mfrow=c(2,2))

#TCGA
TCGA.umap = umap(t(TCGA_genExp))
plot(TCGA.umap$layout,pch=16,xlab = "",ylab = "", main = "TCGA")

TCGA.umap.scale = umap(scale(t(TCGA_genExp)))
plot(TCGA.umap.scale$layout,pch=16,xlab = "",ylab = "", main = "TCGA scaled")


#Lim
Lim.umap = umap(t(Lim_genExp1))
plot(Lim.umap$layout,pch=16,xlab = "",ylab = "", main = "Lim")

Lim.umap.scale = umap(scale(t(Lim_genExp1)))
plot(Lim.umap.scale$layout,pch=16,xlab = "",ylab = "", main = "Lim scaled")


#Dir
Dir.umap = umap(t(Dir_genExp))
plot(Dir.umap$layout,pch=16,xlab = "",ylab = "", main = "Director")

Dir.umap.scale = umap(scale(t(Dir_genExp)))
plot(Dir.umap.scale$layout,pch=16,xlab = "",ylab = "", main = "Director scaled")


#Combine
common_genes<-intersect(rownames(Dir_genExp),rownames(Lim_genExp1))
common_genes<-sort(intersect(common_genes,rownames(TCGA_genExp)))

rbind.data.frame(t(TCGA_genExp[common_genes,]),
                 t(Lim_genExp1[common_genes,]),
                 t(Dir_genExp[common_genes,])) ->genExp

rbind.data.frame(scale(t(TCGA_genExp[common_genes,]),center = T, scale = T),
                 scale(t(Lim_genExp1[common_genes,]),center = T, scale = T),
                 scale(t(Dir_genExp[common_genes,]),center = T, scale = T)) ->genExp_scaled

study<-as.factor(c(rep("TCGA",dim(TCGA_genExp)[2]),
         rep("Lim",dim(Lim_genExp1)[2]),
         rep("Director",dim(Dir_genExp)[2])))

ALL.umap<-umap(genExp)
plot(ALL.umap$layout,col=study,pch=16,xlab = "",ylab = "")

ALL.scaled.umap<-umap(genExp_scaled)
plot(ALL.scaled.umap$layout,col=study,pch=16,xlab = "",ylab = "")

dev.off()

############### 3) Make merged pheno file ##################
"tcga, lim, dir full exp and pheno combined scale only vs combat"
"pheno data sample id, dataset, detailed dataset, os, time, 
age, sex, race, smoking, stage, recurrence/relapse
group, confidence"
"dbu matrix, dbu group assignments preprocessed"
"for stage conversion used ACS table for TNM to AJCC stages
https://www.cancer.org/cancer/lung-cancer/detection-diagnosis-staging/staging-nsclc.html"

#make combined pheno file
TCGA_pheno_merge<-data.frame(Study = "TCGA",
                             Study2 = "TCGA",
                             samp_id = as.vector(TCGA_pheno$sample), 
                             OS_status = ifelse(TCGA_pheno$vital_status == "Dead", 1, 0),
                             OS_months = as.vector(ifelse(TCGA_pheno$vital_status == "Dead",TCGA_pheno$days_to_death,TCGA_pheno$days_to_last_follow_up)/30.4),
                             Age = TCGA_pheno$age_at_index,
                             Sex = TCGA_pheno$gender,
                             Race = recode(TCGA_pheno$race,
                                           `american indian or alaska native` = "Other",
                                           `asian` = "Other",
                                           `black or african american` = "AA",
                                           `not reported` = "NA",
                                           `white` = "White"),
                             Smoker = recode_factor(TCGA_pheno$paper_Smoking.Status,
                                                     `Current reformed smoker for > 15 years` = "Ex",
                                                    `Current reformed smoker for < or = 15 years` = "Ex",
                                                    `Current smoker` = "Current",
                                                    `Lifelong Non-smoker` = "Never",
                                                    `[Not Available]` = "NA"),
                             Stage = recode_factor(TCGA_pheno$ajcc_pathologic_stage,
                                                   `Stage I` = "I",
                                                   `Stage IA` = "I",
                                                   `Stage IB` = "I",
                                                   `Stage II` = "II",
                                                   `Stage IIA` = "II",
                                                   `Stage IIB` = "II",
                                                   `Stage IIIA` = "III",
                                                   `Stage IIIB` = "III",
                                                   `Stage IV` = "IV"))
Lim_pheno_merge<-data.frame(Study = "Lim",
                             Study2 = as.vector(Lim_pheno1$Dataset),
                             samp_id = as.vector(rownames(Lim_pheno1)), 
                             OS_status = Lim_pheno1$Overall.survival..0..alive.1.deceased.,
                             OS_months = as.vector(Lim_pheno1$Overall.survival..month.),
                            Age = Lim_pheno1$Age,
                            Sex = ifelse(as.vector(Lim_pheno1$Gender)=="na", NA, as.vector(Lim_pheno1$Gender)),
                            Race = recode(Lim_pheno1$Race,
                                          `African American` = "AA"),
                            Smoker = recode_factor(Lim_pheno1$Smoking,
                                                   `Ever-smoker` = "Ex",
                                                   `Ex-smoker` = "Ex",
                                                   `Former` = "Ex",
                                                   `Never-smoker` = "Never",
                                                   `Non-smoking` = "Never"),
                            Stage = recode_factor(Lim_pheno1$Stage,
                                                  `1A` = "I",
                                                  `1B` = "I",
                                                  `2A` = "II",
                                                  `2B` = "II",
                                                  ` pT2N0` = "II",
                                                  `1` = "I",
                                                  ` 1A` = "I",
                                                  ` 1B` = "I",
                                                  `2` = "II",
                                                  ` 2A` = "II",
                                                  ` 2B` = "II",
                                                  `3A` = "III",
                                                  `3B` = "III",
                                                  `4` = "IV"))
Dir_pheno_merge<-data.frame(Study = "Director",
                             Study2 = "Director",
                             samp_id = as.vector(Dir_pheno$geo_accession), 
                             OS_status = ifelse(Dir_pheno$`vital_status:ch1` == "Dead", 1, 0),
                             OS_months = Dir_pheno$`months_to_last_contact_or_death:ch1`,
                            Age = Dir_pheno$`age:ch1`,
                            Sex = recode(Dir_pheno$`Sex:ch1`,
                                         Male = "male",
                                         Female = "female"),
                            Race = recode(Dir_pheno$`race:ch1`,
                                          `Asian` = "Other",
                                          `Black or African American` = "AA",
                                          `Native Hawaiian or Other Pacific Islander` = "Other",
                                          `Not Reported` = "NA",
                                          `Unknown` = "NA"),
                            Smoker = recode_factor(Dir_pheno$`smoking_history:ch1`,
                                                   `--` = "NA",
                                                   `Currently smoking` = "Current",
                                                   `Never smoked` = "Never",
                                                   `Smoked in the past` = "Ex",
                                                   `Unknown` = "NA"),
                            Stage = recode_factor(Dir_pheno$`disease_stage:ch1`,
                                                  `pN0pT1` = "I",
                                                  `pN0pT2` = "I/II",
                                                  `pN0pT3` = "II",
                                                  `pN0pT4` = "III",
                                                  `pN1pT1` = "II",
                                                  `pN1pT2` = "II",
                                                  `pN1pT3` = "III",
                                                  `pN1pT4` = "III",
                                                  `pN2pT1` = "III",
                                                  `pN2pT2` = "III",
                                                  `pN2pT3` = "III",
                                                  `pN2pT4` = "III",
                                                  `pNXpT1` = "NA",
                                                  `pp` = "NA"))

pheno_merged<-rbind.data.frame(TCGA_pheno_merge,
                               Lim_pheno_merge,
                               Dir_pheno_merge)
pheno_merged$Smoker<-ifelse(pheno_merged$Smoker == "NA", NA, as.vector(pheno_merged$Smoker))
pheno_merged$Race<-ifelse(pheno_merged$Race == "NA", NA, as.vector(pheno_merged$Race))
pheno_merged$Stage<-ifelse(pheno_merged$Stage == "NA", NA, as.vector(pheno_merged$Stage))
pheno_merged$samp_id<-as.vector(pheno_merged$samp_id)

############## 4) plots ##############

#plot pca of combined data
mypca<-prcomp(genExp_scaled,scale.=FALSE)

library(ggplot2)
ggplot(data = data.frame(mypca$x),
       aes(x=PC1,y=PC2,
           colour=as.factor(pheno_merged$Study)))+
  geom_point(size=0.5)+
  theme_classic()
ggsave("./Results/TCGA_Lim_Dir_pca.png",width = 5,height=3)

############### 4) SAVE ##################
save(genExp_scaled, pheno_merged,file = "./data/Merged_genexp_pheno.Rdata")
save(TCGA_genExp,Lim_genExp1,Dir_genExp,
     TCGA_pheno, Lim_pheno1, Dir_pheno,file = "./data/TCGA_Lim_Dir_genexp_pheno.Rdata")




