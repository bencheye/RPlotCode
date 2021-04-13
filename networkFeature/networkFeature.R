
library(dplyr)
library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)

fileList <- list.files('Input/')

dat1 <- read.table(paste0('Input/', fileList[1]), header = T, sep = '\t')
dat2 <- read.table(paste0('Input/', fileList[2]), header = T, sep = ',')
dat3 <- read.table(paste0('Input/', fileList[3]), header = T, sep = ',')
dat4 <- read.table(paste0('Input/', fileList[4]), header = T, sep = ',')

dat <- dat1
preprocess <- function(dat){
  tmpDat <- data.frame()
  for(i in 1:(dim(dat)[2]-1)){
    tmp <- data.frame(rep(colnames(dat)[i+1], 10), dat$Portion, dat[ ,i+1])
    tmpDat <- rbind(tmpDat, tmp)
  }
  colnames(tmpDat) <- c('Level', 'Portion', 'Value')
  tmpDat$Level <- factor(tmpDat$Level, levels = colnames(dat)[2:dim(dat)[2]])
  tmpDat$Portion <- factor(tmpDat$Portion)
  return(tmpDat)
}

dat1 <- preprocess(dat1)
dat2 <- preprocess(dat2)
dat3 <- preprocess(dat3)
dat4 <- preprocess(dat4)
# theme -------------------------------------------------------------------
mycol<-c("#bcbd22","#4ad693","#1f77b4","#cde64c","#9ac901","#56028e","#a2d1cf","#c8b631","#6dbceb","#c24642")

mytheme<-theme_bw()+theme(
  panel.grid=element_blank(),
  text=element_text(colour="black"),
  axis.text=element_text(size=rel(1)),
  axis.title=element_text(size=rel(1.25),lineheight=0.4,hjust=0.5),
  panel.border=element_rect(size=0.3,colour="black"),
  plot.title=element_text(size=rel(1.5),lineheight=0.7,hjust=0.5),
  #legend.title=element_text(colour=""),
  legend.key=element_rect(colour = "white"),
  axis.line=element_blank()
)

summPlot <- function(dat, numPvalue, legend = 'True'){
  #num <- sum(dat$Value)/length(unique(dat$Level))
  p <- ggplot(data = dat, mapping = aes(x = Level, y  = Value, fill = Portion)) + 
    geom_bar(stat = 'identity', position = 'dodge') + labs(x='', y='') + mytheme +
    scale_fill_manual(breaks = factor(dat$Portion), values=rev(mycol)) + 
    geom_hline(aes(yintercept = numPvalue), colour = 'red', linetype='dashed')
  if(legend == 'Not'){
    p <- p + guides(fill=FALSE)
  }
  return(p)
}

p1 <- summPlot(dat1, 17)
p2 <- summPlot(dat2, 16, legend = 'True')
p3 <- summPlot(dat3, 12.5)
p4 <- summPlot(dat4, 17)

ggpubr::ggarrange(p1, p2, p3, p4, labels = c("A", "B", "C", "D"), 
                  ncol = 2, nrow = 2)
ggsave('ori-driver-networkFeature_result_20210409.pdf', width = 12, height = 6)
ggsave('ori-driver-networkFeature_result_20210409.tiff', width = 12, height = 6)
