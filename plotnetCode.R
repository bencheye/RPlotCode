###LUAD_DEG_neighbor_network_simplify_edge
###time_2019.07.16
options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(ggpubr)

# Input -------------------------------------------------------------------
fileName <- list.files('Input/')
#
oriNetTop <- read.csv(paste0('Input/', fileName[1]))
essGenes <- read.csv(paste0('Input/', fileName[2]), header = F)
tfGenes <- read.csv(paste0('Input/', fileName[3]))
#
driver1 <- read.csv(paste0('Input/', fileName[4]))
driver2 <- read.csv(paste0('Input/', fileName[5]))
driver3 <- read.csv(paste0('Input/', fileName[6]))
driver4 <- read.csv(paste0('Input/', fileName[7]))
#
deg1 <- read.csv(paste0('Input/', fileName[12]), sep = '\t')
deg2 <- read.csv(paste0('Input/', fileName[13]))
deg3 <- read.csv(paste0('Input/', fileName[14]))
deg4 <- read.csv(paste0('Input/', fileName[15]))

colnames(essGenes)[1] <- "GeneSymbol"
colnames(deg1)[1] <- "GeneSymbol"
colnames(deg2)[1] <- "GeneSymbol"
colnames(deg3)[1] <- "GeneSymbol"
colnames(deg4)[1] <- "GeneSymbol"
# normalize cytoscape topology output  ------------------------------------
normal_cyto <- function(cyto_topo){
  cyto_topo <- cyto_topo[ ,c(9, 6, 3:5, 19)]
  colnames(cyto_topo) <- c('GeneSymbol', 'Degree', 'Betweenness', 'Closeness', 
                           'ClusterCoef', 'TopologicalCoef')
  cyto_topo['TF'] <- cyto_topo['Essential'] <- 'NOT'
  cyto_topo$Essential[cyto_topo$GeneSymbol %in% essGenes$GeneSymbol] <- 'Essential'
  cyto_topo$TF[cyto_topo$GeneSymbol %in% tfGenes$GeneSymbol] <- 'TF'
  return(cyto_topo)
}

# preprocess plotData -----------------------------------------------------
preprocess_PlotDat <- function(cyto_topo, driverGene, degList){
  cyto_topo['DEG'] <- cyto_topo['DriverGene']<- 'NOT'
  cyto_topo$DriverGene[cyto_topo$GeneSymbol %in% driverGene$GeneSymbol] <- 'DriverGene'
  cyto_topo$DEG[cyto_topo$GeneSymbol %in% degList$GeneSymbol] <- 'DEG'
  #cyto_topo$NOTDriver[!(cyto_topo$GeneSymbol %in% driverGene$GeneSymbol)] <- 'NOTDriver'
  outDat <- datFormTransform(cyto_topo)
  return(outDat)
}

# kuang data transformed into Long Data -----------------------------------
oriDat <- cyto_topo
datFormTransform <- function(oriDat){
  outDat <- data.frame()
  # choose driver genes' topology centrality
  for(varName in colnames(oriDat)[7:dim(oriDat)[2]]){
    tmpDat <- (oriDat[ ,varName] == varName) %>% oriDat[., 2:6]
    tmpDat['Group'] <- varName
    outDat <- rbind(outDat, tmpDat)
  }
  outDat$Group <- factor(outDat$Group, levels = c('DriverGene', 'Essential', 
                                                     'TF', 'DEG'))
  return(outDat)
}

# plotbox
compareBoxPlot <- function(cyto_topo, centralityName){
  #不同检验比较分组
  my_comparisons <- list(c("DriverGene", "Essential"), c("DriverGene","TF"), 
                         c("DriverGene", "DEG"))
  p <- ggboxplot(cyto_topo, x = "Group", y = centralityName,
            color="Group", palette = "jco") +
    stat_compare_means(comparisons = my_comparisons) +
                       #method.args = list(alternative = "greater")) +
    mytheme + labs(x = '') + theme(legend.position="none")
  return(p)
}
# theme
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

# mainFun -----------------------------------------------------------------
mainFun <- function(cyto_topo, driver, deg){
  # different dataset
  plotDat <- preprocess_PlotDat(cyto_topo, driver, deg)
  p1 <- compareBoxPlot(plotDat, 'Degree')
  p2 <- compareBoxPlot(plotDat, 'Betweenness')
  p3 <- compareBoxPlot(plotDat, 'Closeness')
  p4 <- compareBoxPlot(plotDat, 'ClusterCoef')
  p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
  ggpubr::ggarrange(p1, p2, p3, p4, p5, labels = c("A", "B", "C", "D", "E"), 
                    ncol = 2, nrow = 3)
}

# Preprocess --------------------------------------------------------------
cyto_topo <- normal_cyto(cyto_topo = oriNetTop)
# different dataset
mainFun(cyto_topo, driver1, deg1)
ggsave('Lung-ori-driver-functionGene_compare_20210413.pdf', width = 12, height = 12)
ggsave('Lung-ori-driver-functionGene_compare_20210413.tiff', width = 12, height = 12)
mainFun(cyto_topo, driver2, deg2)
ggsave('Breast-ori-driver-functionGene_compare_20210413.pdf', width = 12, height = 12)
ggsave('Breast-ori-driver-functionGene_compare_20210413.tiff', width = 12, height = 12)
mainFun(cyto_topo, driver3, deg3)
ggsave('Prostate-ori-driver-functionGene_compare_20210413.pdf', width = 12, height = 12)
ggsave('Prostate-ori-driver-functionGene_compare_20210413.tiff', width = 12, height = 12)
mainFun(cyto_topo, driver4, deg4)
ggsave('CRC-ori-driver-functionGene_compare_20210413.pdf', width = 12, height = 12)
ggsave('CRC-ori-driver-functionGene_compare_20210413.tiff', width = 12, height = 12)




