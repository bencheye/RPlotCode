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
#
driver1 <- read.csv(paste0('Input/', fileName[4]))
driver2 <- read.csv(paste0('Input/', fileName[5]))
driver3 <- read.csv(paste0('Input/', fileName[6]))
driver4 <- read.csv(paste0('Input/', fileName[7]))
#
degn1 <- read.csv(paste0('Input/', fileName[8]))
degn2 <- read.csv(paste0('Input/', fileName[9]))
degn3 <- read.csv(paste0('Input/', fileName[10]))
degn4 <- read.csv(paste0('Input/', fileName[11]))

# normalize cytoscape topology output  ------------------------------------
normal_cyto <- function(cyto_topo){
  cyto_topo <- cyto_topo[ ,c(9, 6, 3:5)]
  colnames(cyto_topo) <- c('GeneSymbol', 'Degree', 'Betweenness', 'Closeness', 
                           'ClusterCoef')
  cyto_topo['lung'] <- cyto_topo['breast'] <- cyto_topo['prostate'] <- 
    cyto_topo['CRC'] <- 'Non-Driver'
  cyto_topo$lung[cyto_topo$GeneSymbol %in% driver1$GeneSymbol] <- 'Driver'
  cyto_topo$breast[cyto_topo$GeneSymbol %in% driver2$GeneSymbol] <- 'Driver'
  cyto_topo$prostate[cyto_topo$GeneSymbol %in% driver3$GeneSymbol] <- 'Driver'
  cyto_topo$CRC[cyto_topo$GeneSymbol %in% driver4$GeneSymbol] <- 'Driver'
  #cyto_topo['TF'] <- 'NOT'
  #cyto_topo['Essential'] <- 
  #cyto_topo$Essential[cyto_topo$GeneSymbol %in% essGenes$GeneSymbol] <- 'Essential'
  #cyto_topo$TF[cyto_topo$GeneSymbol %in% tfGenes$GeneSymbol] <- 'TF'
  return(cyto_topo)
}

plotDatFun <- function(cyto_topo, centralityType, CancerType){
    plotDat <- cyto_topo[order(cyto_topo[ ,centralityType], decreasing = T), ]
    plotDat['Order'] <- 1:dim(plotDat)[1]
    xTmp <- plotDat[plotDat[CancerType] == 'Driver', ]  
    pp <- ggplot(plotDat, aes_string('Order', centralityType))+labs(x='') +
      geom_point(size = 2.0, shape = 16, color = 'grey') + mytheme +
      geom_point(data = xTmp, aes_string(x= 'Order', y = centralityType), color = "#FC4E07")
    return(pp)
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
# preprocess plotData -----------------------------------------------------
cyto_topo <- normal_cyto(oriNetTop)
degn1 <- normal_cyto(degn1)
degn2 <- normal_cyto(degn2)
degn3 <- normal_cyto(degn3)
degn4 <- normal_cyto(degn4)

# orign
p1 <- plotDatFun(cyto_topo, 'Degree', 'CRC')
p2 <- plotDatFun(cyto_topo, 'Betweenness', 'CRC')
p3 <- plotDatFun(cyto_topo, 'Closeness', 'CRC')
p4 <- plotDatFun(cyto_topo, 'ClusterCoef', 'CRC')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp4 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp4
ggsave('CRC-oriPPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('CRC-oriPPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(cyto_topo, 'Degree', 'lung')
p2 <- plotDatFun(cyto_topo, 'Betweenness', 'lung')
p3 <- plotDatFun(cyto_topo, 'Closeness', 'lung')
p4 <- plotDatFun(cyto_topo, 'ClusterCoef', 'lung')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp1 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp1
ggsave('lung-oriPPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('lung-oriPPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(cyto_topo, 'Degree', 'breast')
p2 <- plotDatFun(cyto_topo, 'Betweenness', 'breast')
p3 <- plotDatFun(cyto_topo, 'Closeness', 'breast')
p4 <- plotDatFun(cyto_topo, 'ClusterCoef', 'breast')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp2 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp2
ggsave('breast-oriPPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('breast-oriPPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(cyto_topo, 'Degree', 'prostate')
p2 <- plotDatFun(cyto_topo, 'Betweenness', 'prostate')
p3 <- plotDatFun(cyto_topo, 'Closeness', 'prostate')
p4 <- plotDatFun(cyto_topo, 'ClusterCoef', 'prostate')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp3 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp3
ggsave('prostate-oriPPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('prostate-oriPPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

ppOri <- ggpubr::ggarrange(pp1, pp2, pp3, pp4, ncol = 2, nrow = 2,
                           labels = c('A', 'B', 'C', 'D'))
ppOri
ggsave('Summary-oriPPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 20, height = 15)
ggsave('Summary-oriPPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 20, height = 15)


# DEGN PPI  ---------------------------------------------------------------
p1 <- plotDatFun(degn1, 'Degree', 'lung')
p2 <- plotDatFun(degn1, 'Betweenness', 'lung')
p3 <- plotDatFun(degn1, 'Closeness', 'lung')
p4 <- plotDatFun(degn1, 'ClusterCoef', 'lung')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp1 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp1
ggsave('lung-DEGN-PPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('lung-DEGN-PPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(degn2, 'Degree', 'breast')
p2 <- plotDatFun(degn2, 'Betweenness', 'breast')
p3 <- plotDatFun(degn2, 'Closeness', 'breast')
p4 <- plotDatFun(degn2, 'ClusterCoef', 'breast')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp2 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp2
ggsave('breast-DEGN-PPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('breast-DEGN-PPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(degn3, 'Degree', 'prostate')
p2 <- plotDatFun(degn3, 'Betweenness', 'prostate')
p3 <- plotDatFun(degn3, 'Closeness', 'prostate')
p4 <- plotDatFun(degn3, 'ClusterCoef', 'prostate')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp3 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp3
ggsave('prostate-DEGN-PPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('prostate-DEGN-PPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

p1 <- plotDatFun(degn4, 'Degree', 'CRC')
p2 <- plotDatFun(degn4, 'Betweenness', 'CRC')
p3 <- plotDatFun(degn4, 'Closeness', 'CRC')
p4 <- plotDatFun(degn4, 'ClusterCoef', 'CRC')
#p5 <- compareBoxPlot(plotDat, 'TopologicalCoef')
pp4 <- ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
pp4
ggsave('CRC-DEGN-PPI-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 12, height = 10)
ggsave('CRC-DEGN-PPI-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 12, height = 10)

DEGN <- ggpubr::ggarrange(pp1, pp2, pp3, pp4, ncol = 2, nrow = 2,
                           labels = c('A', 'B', 'C', 'D'))
DEGN
ggsave('Summary-DEGN-driverGene-networkTopologyDistribution_20210420.pdf', 
       width = 20, height = 15)
ggsave('Summary-DEGN-driverGene-networkTopologyDistribution_2021042.tiff',
       width = 20, height = 15)
