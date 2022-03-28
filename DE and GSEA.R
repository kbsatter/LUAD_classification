## lung adenocarcinoma differential expressiom
## other characteristics
## immune cell infiltration
## mutation types
## methylation types

library(tidyverse)
library(magrittr)
library(limma)
library(EnhancedVolcano)
library(ComplexHeatmap)
###
pheno = read.csv("TCGA.LUNG.sampleMap_LUNG_clinicalMatrix", sep = '\t')
pheno %<>% filter(X_primary_disease == "lung adenocarcinoma")
pheno$sampleID = gsub("-", ".", pheno$sampleID)
df = read.csv("TCGA.LUAD.sampleMap_HiSeqV2.gz", sep = "\t", row.names = 1)
match(colnames(df), pheno$sampleID)
pheno = pheno[ match(colnames(df), pheno$sampleID), ]
table(pheno$sample_type)
pheno %<>% mutate(stype = ifelse(sample_type == "Solid Tissue Normal",
                                 "STN", "Tumor"))
pheno2 = read.csv("LUAD_Subgroup.csv")
table(pheno2$LUAD_SG)
pheno2 %<>% filter(LUAD_SG != 5)

Subgroup = data.frame(pheno2$sampleID, pheno2$LUAD_SG)
colnames(Subgroup)[1] = "sampleID"
pheno3 = left_join(pheno, Subgroup, by.x = sampleID)

#### limma
pheno3$pheno2.LUAD_SG %<>% replace_na("STN")
ctype = pheno3$pheno2.LUAD_SG
lung = cbind.data.frame(ctype, t(df) %>% as.data.frame())

design = model.matrix( ~ 0 + lung$ctype)
colnames(design) = c("TP1", "TP2", "TP3", "TP4", "STN")
fit = lmFit(t(lung[ , -1]), design)
contrast.matrix = makeContrasts(
  "TP1" = TP1 - STN,
  "TP2" = TP2 - STN,
  "TP3" = TP3 - STN,
  "TP4" = TP4 - STN,
  levels = design
)
fit2 = contrasts.fit(fit, contrast.matrix)
fit2 = eBayes(fit2)
res.TP1 = topTable(fit2, coef = "TP1", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP2 = topTable(fit2, coef = "TP2", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP3 = topTable(fit2, coef = "TP3", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP4 = topTable(fit2, coef = "TP4", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")

####
## gsea
library(org.Hs.eg.db)
library(fgsea)
pathways = gmtPathways("c2.cp.reactome.v7.5.1.symbols.gmt")
# ranks = res.TP1$logFC
# names(ranks) = rownames(res.TP1)
# TP1.GS = fgsea(pathways, stat = ranks,
#       minSize = 15, maxSize = 500)
# TP1.GS %<>% filter(padj <= 0.05)

franks = function(x){
  x %<>% arrange(-logFC)
  ranks = x$t
  names(ranks) = rownames(x)
  res.fgsea = fgsea(pathways,
                    stat = ranks,
                    minSize = 15,
                    maxSize = 500,
                    eps = 0)
  res.fgsea %<>% filter(padj <= 0.05) %>% arrange(-NES)
}

TP1.GSEA = franks(res.TP1)
TP2.GSEA = franks(res.TP2)
TP3.GSEA = franks(res.TP3)
TP4.GSEA = franks(res.TP4)
TP1.GSEA = apply(TP1.GSEA, 2, as.character)
write.csv(TP1.GSEA, "TP1 Vs N.csv")
write.csv(apply(TP2.GSEA, 2, as.character), "TP2 Vs N.csv")
write.csv(apply(TP3.GSEA, 2, as.character), "TP3 Vs N.csv")
write.csv(apply(TP4.GSEA, 2, as.character), "TP4 Vs N.csv")

### difference between the types

contrast.matrix2 = makeContrasts(
  "TP1" = TP1 - (TP2 + TP3 + TP4)/3,
  "TP2" = TP2 - (TP1 + TP3 + TP4)/3,
  "TP3" = TP3 - (TP1 + TP2 + TP4)/3,
  "TP4" = TP4 - (TP1 + TP2 + TP3)/3,
  levels = design
)
fit3 = contrasts.fit(fit, contrast.matrix2)
fit3 = eBayes(fit3)
res.TP1.1 = topTable(fit3, coef = "TP1", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP2.1 = topTable(fit3, coef = "TP2", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP3.1 = topTable(fit3, coef = "TP3", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")
res.TP4.1 = topTable(fit3, coef = "TP4", n = dim(lung)[2]- 1,
                   p.value = 0.05, adjust.method = "fdr")

TP1.GSEA.1 = franks(res.TP1.1)
TP2.GSEA.1 = franks(res.TP2.1)
TP3.GSEA.1 = franks(res.TP3.1)
TP4.GSEA.1 = franks(res.TP4.1)

#### write the files
write.csv(apply(TP1.GSEA.1, 2, as.character), "TP1 Vs other.csv")
write.csv(apply(TP2.GSEA.1, 2, as.character), "TP2 Vs other.csv")
write.csv(apply(TP3.GSEA.1, 2, as.character), "TP3 Vs other.csv")
write.csv(apply(TP4.GSEA.1, 2, as.character), "TP4 Vs other.csv")



###
mcth1 = TP2.GSEA.1[TP2.GSEA.1$pathway %in% TP1.GSEA.1$pathway, ]
mcth2 = mcth1[mcth1$pathway %in% TP3.GSEA.1$pathway, ]
mcth3 = mcth2[mcth2$pathway %in% TP4.GSEA.1$pathway, ]

###
TP1.GSEA.1 %<>% add_column(Type = rep("TP1", dim(TP1.GSEA.1)[1]))
TP2.GSEA.1 %<>% add_column(Type = rep("TP2", dim(TP2.GSEA.1)[1]))
TP3.GSEA.1 %<>% add_column(Type = rep("TP3", dim(TP3.GSEA.1)[1]))
TP4.GSEA.1 %<>% add_column(Type = rep("TP4", dim(TP4.GSEA.1)[1]))

TP1.GSEA.1.2 = TP1.GSEA.1[match(mcth3$pathway, TP1.GSEA.1$pathway), ]
TP2.GSEA.1.2 = TP2.GSEA.1[match(mcth3$pathway, TP2.GSEA.1$pathway), ]
TP3.GSEA.1.2 = TP3.GSEA.1[match(mcth3$pathway, TP3.GSEA.1$pathway), ]
TP4.GSEA.1.2 = TP4.GSEA.1[match(mcth3$pathway, TP4.GSEA.1$pathway), ]

Gdata = rbind.data.frame(TP1.GSEA.1.2, TP2.GSEA.1.2, TP3.GSEA.1.2, TP4.GSEA.1.2)

Gdata$pathway = str_remove(Gdata$pathway, "REACTOME_")

Gdata %>% 
  ggplot(aes(x = Type, y = pathway, size = NES, color = log2(padj))) +
  geom_point() + theme_classic() + theme(text = element_text(size = 14, color = "black")) +
  scale_color_gradient(low = "red", high = "blue", space = "Lab",
                       aesthetics = "color")
