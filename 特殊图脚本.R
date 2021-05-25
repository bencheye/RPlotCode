### 画出KEGG通路富集结果

#特定通路作图
# pathway
# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=gene_list, pathway.id="hsa05012", species = 'hsa')
# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=gene_list, pathway.id="hsa05012", species = 'hsa', kegg.native = F)
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
                     knitr::include_graphics("hsa05012.pathview.png")

### 画PuMed研究趋势图
# pubmed趋势图
library(europepmc)
terms <- kk2$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)
