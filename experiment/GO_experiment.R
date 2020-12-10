rm(list=ls())
# from http://yulab-smu.top/clusterProfiler-book/chapter5.html
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
head(gene.df)
dim(gene.df)

################

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

head(ggo)

################

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

?clusterProfiler::enrichGO

#################################

library(org.Mmu.eg.db)
xx <- as.list(org.Mmu.egALIAS2EG)
xx <- xx[!is.na(xx)]
# oh.. this is everything I need I think.
gene_vec <- read.csv("../data-raw/All_human_genes.txt", header = F)
entrez_id <- sapply(1:nrow(gene_vec), function(i){
  if(i %% floor(nrow(gene_vec)/10) == 0) cat('*')
  
  id <- gene_vec[i,1]
  idx <- which(names(xx) == id)
  
  if(length(idx) == 0){
    return(NA)
  } else {
    xx[[idx[1]]]
  }
})
sum(is.na(entrez_id))
entrez_vec <- unlist(entrez_id)
entrez_vec <- entrez_vec[!is.na(entrez_vec)]

ego2 <- enrichGO(gene         = entrez_vec[1:50],
                universe      = entrez_vec,
                OrgDb         = org.Mmu.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
ego2

