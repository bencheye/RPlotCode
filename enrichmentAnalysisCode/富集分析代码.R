###Driver genes enrichment analysis
###
###
#BiocManager::install(c("gridExtra"))
#install.packages('org.Hs.eg.db')
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(dplyr)
library(gridExtra)
library(ggplot2)
options(stringsAsFactors = F)

path <- "/Users/benche/analysis/MyPaperFigure/enrichment"
setwd(path)

# Input -------------------------------------------------------------------
driverList <- list.files('Input/')
driver1 <- read.csv(paste0('Input/', driverList[1]))
driver2 <- read.csv(paste0('Input/', driverList[2]))
driver3 <- read.csv(paste0('Input/', driverList[3]))
driver4 <- read.csv(paste0('Input/', driverList[4]))

# transform GeneSymbol into GeneID
eg1 <- bitr(driver1[ ,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg2 <- bitr(driver2[ ,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg3 <- bitr(driver3[ ,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
eg4 <- bitr(driver4[ ,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# theme -------------------------------------------------------------------
mytheme<-theme_bw()+theme(
  panel.grid=element_blank(),
  text=element_text(colour="black"),
  axis.text=element_text(size=rel(1)),
  axis.title=element_text(size=rel(1.25),lineheight=0.4,hjust=0.5),
  panel.border=element_rect(size=0.3,colour="black"),
  plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
  legend.title=element_text(colour="white"),
  legend.key=element_rect(colour = "white"),
  axis.line=element_blank()
)


# BP ----------------------------------------------------------------------
eg <- eg1
BpPlot <- function(eg){
  ego_bp <- enrichGO(gene          = eg$ENTREZID,
                      #universe      = eg$SYMBOL,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2,
                      readable      = TRUE)
  res <- data.frame(ego_bp@result)
  #res <- data.frame(ego_bp1@result)
  numA <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][1])}) %>% as.numeric()
  numB <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][2])}) %>% as.numeric()
  res$GeneRatio <- round(numA/numB, 2)
  res <- res[order(res$Count[1:10]), ]
  res$Description <- factor(res$Description, levels = res$Description)
  p <- ggplot(data = res[1:10, ], 
         mapping = aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) + 
    geom_point() + mytheme + labs(y='') + 
    scale_color_gradient(low="#1f77b4", high="#cde64c")
  return(p)
}
p1 <- BpPlot(eg1)
p2 <- BpPlot(eg2)
p3 <- BpPlot(eg3)
p4 <- BpPlot(eg4)
#titles <- c('Lung Cancer', 'Breast Cancer', 'Prostate Cancer', 'CRC')
ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), 
                  ncol = 2, nrow = 2)
ggsave('GO_enrichment_result_20210409.pdf', width = 15, height = 8)
ggsave('GO_enrichment_result_20210409.tiff', width = 15, height = 8)
# KEGG --------------------------------------------------------------------
KeggPlot <- function(eg){
  kk_bp <- enrichKEGG(eg$ENTREZID, organism="hsa",
                   keyType = "ncbi-geneid",
                   pvalueCutoff=0.05, pAdjustMethod="BH",
                   qvalueCutoff=0.1)
  res <- data.frame(kk_bp@result)
  #res <- data.frame(ego_bp1@result)
  numA <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][1])}) %>% as.numeric()
  numB <- sapply(res$GeneRatio, function(x){return(strsplit(x, '/')[[1]][2])}) %>% as.numeric()
  res$GeneRatio <- round(numA/numB, 2)
  res <- res[order(res$Count[1:10]), ]
  res$Description <- factor(res$Description, levels = res$Description) 
  p <- ggplot(data = res[1:10, ], 
              mapping = aes(x = GeneRatio, y = Description, colour = p.adjust, size = Count)) + 
    geom_point() + mytheme + labs(y='') +
    scale_color_gradient(low="#1f77b4", high="#cde64c")
  return(p)
}
pp1 <- KeggPlot(eg1)
pp2 <- KeggPlot(eg2)
pp3 <- KeggPlot(eg3)
pp4 <- KeggPlot(eg4)
ggpubr::ggarrange(pp1, pp2, pp3, pp4, labels = c("A", "B", "C", "D"), 
                  ncol = 2, nrow = 2)
ggsave('KEGG_enrichment_result_20210409.pdf', width = 15, height = 8)
ggsave('KEGG_enrichment_result_20210409.tiff', width = 15, height = 8)

