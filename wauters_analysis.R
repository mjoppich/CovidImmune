library("Seurat")
library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)
library(tidyverse)
library(pbapply)
library(clusterProfiler)
library("enrichplot")
library(data.table)
library("writexl")
library("EnhancedVolcano")
library(ggsunburst)
library(ggrepel)
library(RColorBrewer)
library(svglite)
library(ggpubr)
library("glmGamPoi")

save_plot = function(plotobj,outname, fig.width, fig.height)
{
  print(paste(outname, fig.width, fig.height))
  
  fname=paste(outname, "png", sep=".")
  print(paste("Saving to file", fname))
  png(filename=fname, width = fig.width, height = fig.height, units = 'in', res = 300)#width = fig.width*100, height=fig.height*100)
  plot(plotobj)
  dev.off()
  
  fname=paste(outname, "pdf", sep=".")
  print(paste("Saving to file", fname))
  pdf(file=fname, width = fig.width, height=fig.height)
  plot(plotobj)
  dev.off()
  

  fname=paste(outname, "svg", sep=".")
  print(paste("Saving to file", fname))
  svglite::svglite(file = fname, width = fig.width, height = fig.height)
  plot(plotobj)
  dev.off()
  
  return(plotobj)
}

"""
The count and annotation data were downloaded from the Lambrecht lab website (http://covid19.lambrechtslab.org/) [cite https://www.nature.com/articles/s41422-020-00455-9].

The data were splitted by patient/sample and processed following the SCTransfrom data integration vignette, using 2000 features for data integration and regressing out ribosomal and mitochonrial content, UMI count and S/G2M-Score.
The in the original paper mentioned patient-specific cells (mostly patient BAL019) were identified and corresponding cell IDs were saved.

In a second pass, the previously identified cells were removed. Each patient sample wes normalized, identified 3000 variable genes and scaled, while regressing out ribosomal and mitochonrial content, UMI count and S/G2M-Score.
In order to perform the integration with 3000 features, the reciprocal PCA integration vignette was followed.

Celltype and disease-state annotations were mapped by cell-ID from the original author's annotation.
The analysis was performed with R 4.0.1 using Seurat 4.0.0 [cite].
The statistical comparison between the disease stages was performed using the stat_compare_means function of the ggpubr package ( https://github.com/kassambara/ggpubr ).

"""

toCellPrefix = function(x) {
    datasetprefix = x
    datasetprefix = str_replace_all(datasetprefix, "[./]", "")
    datasetprefix = str_split(datasetprefix, "_")[[1]][1]
    
    return(paste(datasetprefix, sep=""))
}

makeSeuratObj = function(matrix, proj, minUMIs, plots)
{
    obj = CreateSeuratObject(matrix, project=proj)
    print("Renaming Cells")
    obj <- RenameCells(obj, add.cell.id=proj)
    
    print(paste("Seurat obj project", obj@project.name))
    
    #obj[["percent.mtrp"]] <- PercentageFeatureSet(obj, pattern = "^mt-|^Rps|^Rpl")
    
    mtPattern = "^MT-"
    rplPattern = "^RPL"
    rpsPattern = "^RPS"
    rpPattern = "^RPS|^RPL"
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=mtPattern)]
    print(paste("Got a total of mt-Genes:", length(selGenes), paste0(head(unlist(selGenes)), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rplPattern)]
    print(paste("Got a total of Rpl-Genes:", length(selGenes), paste0(head(unlist(selGenes)), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpsPattern)]
    print(paste("Got a total of Rps-Genes:", length(selGenes), paste0(head(unlist(selGenes)), collapse=", ")))
    
    selGenes = rownames(obj)[grepl(rownames(obj), pattern=rpPattern)]
    print(paste("Got a total of Rp-Genes:", length(selGenes), paste0(head(unlist(selGenes)), collapse=", ")))
    
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mtPattern)
    obj[["percent.rpl"]] <- PercentageFeatureSet(obj, pattern = rplPattern)
    obj[["percent.rps"]] <- PercentageFeatureSet(obj, pattern = rpsPattern)
    obj[["percent.rp"]] <- PercentageFeatureSet(obj, pattern = rpPattern)
    
    if (plots)
    {
      plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
      plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      show(plot1 + scale_x_continuous(n.breaks = 20) + scale_y_continuous(n.breaks = 20))
    }
    
    # mt content: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6072887/
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > minUMIs & percent.mt < 15 & percent.rp < 40)
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE, nfeatures = 3000)
    
    return(obj)
}



iaBalf = readRDS("../data2/Allcells.counts.rds")
iaBalfAnnot = read.csv("../data2/Allcells.meta.data.csv")
removeCells =read.csv("remove_cells_bal19")

allSamples = unique(unlist(lapply(str_split(colnames(iaBalf), "_"), function(x){return(x[1])})))
objlist = list()

for (sample in allSamples)
{

    cellnames = grep(x=colnames(iaBalf), pattern=sample, value=T)
    sampleMatrix = iaBalf[, cellnames]

    newColNames = unlist(lapply(str_split(colnames(sampleMatrix), "_"), function(x){return(x[2])}))
    print(head(newColNames))
    colnames(sampleMatrix) = newColNames

    acceptCellnames = setdiff(colnames(sampleMatrix), removeCells$cells)
    sampleMatrix = sampleMatrix[, acceptCellnames]

    expName = paste("iabalf_", sample, sep="")
    print(paste(sample, expName, length(cellnames), dim(sampleMatrix)))
    sampleObj = makeSeuratObj(sampleMatrix, expName, 300, F)

    print(head(colnames(sampleObj)))

    objlist[[sample]] = sampleObj
}


print("cells per experiment")
print(mapply(sum, lapply(objlist, function(x) {dim(x)[2]})))
print("total cells")
print(sum(mapply(sum, lapply(objlist, function(x) {dim(x)[2]}))))

objlist <- lapply(X = objlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = objlist, nfeatures = 3000)
objlist <- lapply(X = objlist, FUN = function(x) {

    print(x)

      s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes

    x <- CellCycleScoring(
      x,
      g2m.features = g2m.genes,
      s.features = s.genes)


    x <- ScaleData(x, features = features, verbose = FALSE,vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score"))
    x <- RunPCA(x, features = features, verbose = FALSE)
})

#objlist <- lapply(X = objlist, FUN = function(x) {
#    print(x)
    #x <- ScaleData(x, features = noRPMTFeatures, verbose = FALSE, vars.to.regress = c('percent.rp', 'percent.mt'))
    #x <- RunPCA(x, features = noRPMTFeatures, npcs = 50, verbose = FALSE)
    #suppressWarnings( x <- SCTransform(x, vars.to.regress = c('percent.rp', 'percent.mt', "nCount_RNA","S.Score", "G2M.Score"), method = "glmGamPoi", verbose=FALSE) )
#})

objlist.anchors <- FindIntegrationAnchors(object.list = objlist, reduction = "rpca", dims = 1:50,anchor.features = 3000)

save.image(file="before_integration2000.v4.final.RData", compress = F)

obj.integrated <- IntegrateData(anchorset = objlist.anchors, dims = 1:50, verbose=T)

write_rds(obj.integrated, "integrated_object2000.rds", compress="none")
save.image("integrated2000.v4.final.RData", compress = FALSE)


obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:50, return.model = TRUE)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:50)

DefaultAssay(obj.integrated) = "integrated"
obj.integrated <- FindClusters(obj.integrated, resolution = 0.5)
obj.integrated$idents = Idents(obj.integrated)


iaBalfAnnot$XN = paste("iabalf_", iaBalfAnnot$Cell, sep="")
rownames(iaBalfAnnot) = iaBalfAnnot$XN
cellTypeNames = iaBalfAnnot[colnames(obj.integrated),]$CellType
obj.integrated$orig_celltype = iaBalfAnnot[colnames(obj.integrated),]$CellType

p=DimPlot(obj.integrated, pt.size = 0.001, label=T)
save_plot(p, "dimplot2000_umap", fig.width=12, fig.height=8)

write_rds(obj.integrated, "integrated_object2000.rds", compress="none")
save.image("integrated2000.v4.final.RData", compress = FALSE)


cellList = colnames(obj.integrated)
featVec <- vector(mode="character", length=length(cellList))

featVec[grep(x = cellList, pattern = "^iabalf_BAL001")] = "pneu_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL002")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL003")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL009")] = "covid_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL010")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL011")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL012")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL013")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL014")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL015")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL016")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL017")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL018")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL019")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL020")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL021")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL022")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL023")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL024")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL025")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL026")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL027")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL028")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL029")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL030")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL031")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL032")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL033")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL034")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL035")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL036")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL037")] = "covid_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL038")] = "pneu_mild"
featVec[grep(x = cellList, pattern = "^iabalf_BAL039")] = "covid_severe"
featVec[grep(x = cellList, pattern = "^iabalf_BAL040")] = "covid_severe"

obj.integrated$pneu_state=featVec
unique(obj.integrated$pneu_state)

DefaultAssay(obj.integrated) = "RNA"

p=DimPlot(obj.integrated, pt.size = 0.001, label=T)
save_plot(p, "dimplot2000_umap", fig.width=12, fig.height=8)

p=DimPlot(obj.integrated, split.by="orig.ident", ncol=2, pt.size = 0.001, label=T)
save_plot(p, "dimplot2000_umap_ident", fig.width=12, fig.height=36)

p=DimPlot(obj.integrated, split.by="pneu_state", ncol=4, pt.size = 0.1, label=T)
save_plot(p, "dimplot2000_umap_pneu", fig.width=30, fig.height=8)

p=DimPlot(obj.integrated, group.by="orig_celltype", pt.size = 0.001, label=T)
save_plot(p, "dimplot2000_umap_ct", fig.width=18, fig.height=12)

p=DimPlot(obj.integrated, split.by="orig_celltype", ncol=3, pt.size = 0.001, label=T)
save_plot(p, "dimplot2000_umap_orig_celltype", fig.width=18, fig.height=12)

write_rds(obj.integrated, "integrated_object_annot2000.rds", compress="none")
save.image("integrated2000.v4.final.RData", compress = FALSE)


splitFeaturePlot = function(obj, feature, split.by, title, filename)
{
  pds = FeaturePlot(obj, features = feature, reduction = "umap", split.by=split.by, combine=F,min.cutoff=0, max.cutoff=7,order=T)
  pds[[1]] = pds[[1]] + ggtitle(NULL)

  print(paste(length(pds)))
  pdsRange = c(1:(length(pds)-1))

  for (i in pdsRange)
  {
    print(i)
    pds[[i]] = pds[[i]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")  
  }

  #pds[[1]] = pds[[1]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE) + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  #pds[[2]] = pds[[2]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  #pds[[3]] = pds[[3]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue", guide=FALSE)+ theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")
  pds[[length(pds)]] = pds[[length(pds)]] + scale_color_gradient(limits = c(0,7), low = "lightgrey", high = "blue"             ) + theme(axis.text.x = element_text(face="bold", color="#000000", size=14, angle=0), axis.text.y = element_text(face="bold", color="#000000", size=14, angle=0))+labs(x="", y="")

  prow = cowplot::plot_grid(plotlist=pds, label_x = "a", ncol=4, align="hv")
  prow
  # now add the title
  title <- ggdraw() + 
    draw_label(
      title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )

  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  fplot = plot_grid(
    title, prow,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )

  save_plot(fplot, filename, fig.width=25, fig.height=6)

}
splitFeaturePlot(obj.integrated, "CXCL8", title="CXCL8 expression per disease state", split.by="pneu_state", filename="obj2000_feature_cxcl8_gp_pneustate")


p=DimPlot(obj.integrated, group.by="orig_celltype", ncol=1, pt.size = 0.1, label=T, label.size=6)
save_plot(p, "obj2000_celltypes", fig.width=12, fig.height=8)


library(ggsignif)
allStates = as.factor(setdiff(unique(featVec), c("pneu_severe")))
my_comparisons <- combn(levels(allStates),2, simplify = F)
p=VlnPlot(subset(obj.integrated, orig_celltype=="Neutrophil"), features = c("CXCL8"), group.by="pneu_state", ncol = 1, pt.size=0,y.max = 8, cols=c("covid_mild"="#00E5FF", "covid_severe"="#003280", "pneu_mild"="#68C458", "pneu_severe"="#4D9341"))+ylim(c(0,10))+stat_compare_means(comparisons = my_comparisons,label.y = c(7,8,9), method="t.test")
p=p+ggtitle("CXCL8 Expression in Wauters et al. (Neutrophils, t-test)")
save_plot(p, "obj2000_vlnbxp_cxcl8_gp_pneustate", fig.width=7, fig.height=5)

table(subset(obj.integrated, orig_celltype=="Neutrophil")$pneu_state)


p=VlnPlot(obj.integrated, features = c("CXCL8"), group.by="orig_celltype", split.by="pneu_state", ncol = 1, pt.size=0, cols=c("covid_mild"="#00E5FF", "covid_severe"="#003280", "pneu_mild"="#68C458", "pneu_severe"="#4D9341"))
save_plot(p, "obj2000_vln_cxcl8_gp_pneustate", fig.width=12, fig.height=3)


cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

make_descr_label = function(plot, descr)
{
  descrLabel <- ggdraw() + draw_label(descr, fontface='bold', angle = 0)
  
  pe = cowplot::plot_grid(descrLabel, plot, ncol=1, nrow=2, labels=NULL,rel_heights = c(0.1, 1),
                          align = "h", axis = "l")
  
  return(pe)
}

makeSideBySideDotPlot = function(scobj, plotElems, featureGenes = c(""), group.by="cellnames_manual", col.min = -3, col.max = 3, cols = c("grey", "blue"), title="")
{

  plot_list = list()
  plot_orig_list = list()
  allFoundFeatures = c()
  
  ctFractions = list()
  
  scaled=T
  
  elemCount = 0
  for (plotName in names(plotElems))
  {
    print(plotName)
    plotData = plotElems[[plotName]]
    plotCells = plotData$cells
    plotDescr = plotData$label
    
    scobj_subset = subset(scobj, cells=plotCells)
    plotElem_orig = DotPlot(scobj_subset, features=featureGenes, group.by = group.by, col.min = col.min, col.max=col.max, cols = cols)
    
    scTable = table(scobj_subset[[group.by]])
    scDf = as.data.frame(scTable)
    
    scDf$perc = scDf$Freq / sum(scDf$Freq)
    ctFractions[[plotName]] = scDf    
    
    scobj_subset = NULL
    
    allFoundFeatures = unique(c(allFoundFeatures, as.character(plotElem_orig$data$features.plot)))
    
    plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
    plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none")
    
    plot_orig_list[[plotName]] = plotElem_orig
  }
  
  
  if (scaled)
  {
    
    # initialize Values
    for (plotName in names(plot_orig_list))
    {
        plot_orig_list[[plotName]]$data$avg.exp.scaled2 = 0
        
        plot_orig_list[[plotName]]$data$id = as.character(plot_orig_list[[plotName]]$data$id)
        plot_orig_list[[plotName]]$data$features.plot = as.character(plot_orig_list[[plotName]]$data$features.plot)
        
    }
    
    allIDs = c()
    allFeatures = c()
    for (plotName in names(plot_orig_list))
    {
        allIDs = c(allIDs, plot_orig_list[[plotName]]$data$id)
        allFeatures = c(allFeatures, plot_orig_list[[plotName]]$data$features.plot)
    }
    allIDs = unique(allIDs)
    allFeatures = unique(allFeatures)
    
    print(allIDs)
    print(allFeatures)
    
    
    # calculate avg.exp.scaled2 for each feature
    for (featureName in allFoundFeatures)
    {
     
      allUnscaledValues = NULL
      for (plotName in names(plot_orig_list))
      {
        pData = plot_orig_list[[plotName]]$data
        

        missingCelltypes = setdiff(allIDs, unique(pData$id))
      
        
        for (celltype in missingCelltypes)
        {
          for (feature in allFeatures)
          {
            dfrow = data.frame(avg.exp=0.0,pct.exp=0.0,features.plot=feature, id=celltype, avg.exp.scaled=0.0, avg.exp.scaled2=0.0)
            pData = rbind.data.frame(pData, dfrow, stringsAsFactors = F)
            
          }
        }
        
        pData$id = factor(pData$id, levels = allIDs)
        pData$features.plot = factor(pData$features.plot, levels=allFeatures)
        

        plot_orig_list[[plotName]]$data = pData
        
        
        if (is.null(allUnscaledValues))
        {
          allUnscaledValues = data.frame(plotName=pData[ pData$features.plot==featureName, ]$avg.exp)
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
          allUnscaledValues[["plotName"]] = NULL
          
        } else {
          allUnscaledValues[[plotName]] = pData[ pData$features.plot==featureName, ]$avg.exp  
        }
        
      }
      allUnscaledValues$rnames = as.numeric(rownames(allUnscaledValues))
      
      allUnscaledLong = allUnscaledValues %>% gather(Type, Value, names(plot_orig_list))
      allUnscaledLong$Value = scale(allUnscaledLong$Value)
      
      allScaledValues = allUnscaledLong %>% spread(Type, Value) %>% arrange( order(rnames))
      allScaledValues
      
      for (plotName in names(plot_orig_list))
      {

        plotElem_orig = plot_orig_list[[plotName]]
        pData = plotElem_orig$data
  
        
        # https://github.com/satijalab/seurat/issues/2798
        plotElem_orig$layers[[1]] <- NULL # remove original geom_point layer where the color scale is hard-coded to use scaled average expression
        
        origData = plotElem_orig$data[plotElem_orig$data$features.plot==featureName, ]

        plotElem_orig$data[plotElem_orig$data$features.plot==featureName, "avg.exp.scaled2"] = allScaledValues[,plotName]
        plot_orig_list[[plotName]] = plotElem_orig
      }
    }
    
    for (plotName in names(plot_orig_list))
    {
      plotElem_orig = plot_orig_list[[plotName]]
      pData = plotElem_orig$data
      
      pData2 = merge(x=pData,y=ctFractions[[plotName]],by.x="id", by.y="Var1",all.x=TRUE)
      
      plotElem_orig$data <-pData2 %>% mutate(featuren=as.numeric(features.plot), idn=as.numeric(id), percn=log(perc))
      
      print(plotElem_orig$data)
      
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2>col.max, 'avg.exp.scaled2'] = col.max
      plotElem_orig$data[plotElem_orig$data$avg.exp.scaled2<col.min, 'avg.exp.scaled2'] = col.min
      
      pData = plotElem_orig$data
      
      plotElem_orig <- ggplot(pData,aes(x=featuren, y=idn, colour = avg.exp.scaled2, size = pct.exp)) +
        scale_x_continuous(breaks=plotElem_orig$data$featuren, labels=plotElem_orig$data$features.plot) +
        scale_y_continuous(breaks=plotElem_orig$data$idn, labels=plotElem_orig$data$id)+
        geom_rect(mapping=aes(xmin=featuren-.5, xmax=featuren+.5, ymin = idn-0.5, ymax = idn+0.5, fill = percn), alpha = 0.4, linetype="blank") +
        scale_fill_distiller(palette='Spectral')+
        scale_size_continuous(range = c(0, 10))+
        geom_point() +
        scale_color_gradient2(limits=c(col.min, col.max), low = cols[1], mid=cols[2], high = cols[3])+
        guides(color= guide_colourbar(title="Avg. Expression (scaled)"), size=guide_legend(title="Percent Expressing"), fill=guide_colourbar(title="Cell Abundance (log)"))
      
      #plotElem_orig = plotElem_orig + scale_color_gradient(limits=c(col.min, col.max), low = cols[1], high = cols[2])
      plotElem_orig = plotElem_orig + theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line())
              #scale_color_gradient()+

      plot_orig_list[[plotName]] = plotElem_orig
    }
  }
  
  
  
  elemCount = 0
  plot_list=list()
  for (plotName in names(plot_orig_list))
  {
    
    plotElem_orig = plot_orig_list[[plotName]]
    
    elemCount = elemCount + 1
    
    plotElem = plotElem_orig
    
    #+  theme(
              #axis.text.y = element_blank(),
    #          panel.grid.major.x = element_line( size=.1, color="black" ),
              #axis.line.y = element_line(size = 0),
              #axis.ticks.length.y = unit(0, "points")
    #)
  
    pe = make_descr_label(plotElem, plotElems[[plotName]]$label)
  
    plot_list[[plotName]] = pe
  
  }
  
  legendDescr = 'Average Expression'
  if (scaled)
  {
    legendDescr = 'Average Scaled Expression'  
  }
  
  # extract a legend that is laid out horizontally #
  legend_b <- get_legend(
    plotElem_orig + 
      guides(color = guide_colorbar(title = legendDescr, direction="horizontal"))+ #guide_legend(nrow = 1, override.aes = list(size=6)))+
      theme(legend.position = "bottom")
  )
  
  title <- ggdraw() + draw_label(title, fontface='bold')
  
  ap=cowplot::plot_grid(
    plotlist = plot_list,
    labels = NULL,
    nrow=1,
    align = "v", axis="bt"
  )

  fp = cowplot::plot_grid(title, ap, legend_b, ncol = 1, rel_heights = c(.05, 1, 0.1) )
  return(fp)
}

log_both <- function(x){ifelse(x == 0, 0, log(abs(x)) * sign(x))}

log_both_trans <- 
  function(){
    trans_new(name = 'log_both', 
              transform = log_both,
              inverse = log_both) #not clear what `inverse` does
  }


plotElems = list(
  "COVID-mild"=list(cells=names(obj.integrated$pneu_state[obj.integrated$pneu_state == "covid_mild"]), label="COVID-mild"),
  "COVID-severe"=list(cells=names(obj.integrated$pneu_state[obj.integrated$pneu_state == "covid_severe"]), label="COVID-severe"),
  "PNEU-mild"=list(cells=names(obj.integrated$pneu_state[obj.integrated$pneu_state == "pneu_mild"]), label="PNEU-mild"),
  "PNEU-severe"=list(cells=names(obj.integrated$pneu_state[obj.integrated$pneu_state == "pneu_severe"]), label="PNEU-severe")
)

p=makeSideBySideDotPlot(obj.integrated, plotElems, featureGenes = c("CXCL8", "IL6", "IL1A", "IL1B"), group.by = "orig_celltype", col.min = -4, col.max = 4, title="Scaled Average Expression per Disease Stage of selected genes", cols = c("blue","yellow", "red"))
cowplot::save_plot("obj2000_dotplot_ilgenes.pdf", p, base_width = 12, base_height=6)
cowplot::save_plot("obj2000_dotplot_ilgenes.png", p, base_width = 12, base_height=6)
cowplot::save_plot("obj2000_dotplot_ilgenes.svg", p, base_width = 12, base_height=6)


getCellCountDF = function(scobj, prefix="", group_by="orig.ident", select_by="idents", split_by=NULL, relative=F, outname=NULL,show.percent=F)
{
  
  allClusters = as.character(sort(unique(scobj[[select_by]][,])))
  
  if (!is.null(split_by))
  {
    allSplits = as.character(sort(unique(scobj[[split_by]][,])))
  }
  
  cellCounts = list()
  print(allClusters)
  for (clusterID in allClusters)
  {
    print(clusterID)
    
    if (!is.null(split_by))
    {
      clusterList = list()
      clusterList[["cluster"]] = clusterID
      cs = scobj[,scobj[[select_by]] == clusterID]
      clusterList[["all"]] = nrow(cs[[group_by]])
      for (split_cat in allSplits)
      {
          print(paste(clusterID, split_cat))
        
          cs = scobj[,scobj[[select_by]] == clusterID & scobj[[split_by]] == split_cat]
          allElems = table(cs[[group_by]])
          cs = NULL
        
          for (grp in names(allElems))
          {
            if (!is.null(prefix))
            {
              ngrp = paste(prefix, paste(split_cat, grp, sep="."), sep=".")
  
            } else {
              ngrp = paste(split_cat, grp, sep=".")
            }
            clusterList[[ngrp]] = allElems[[grp]]
          }
      }
      
    } else {
      
        cs = scobj[,scobj[[select_by]] == clusterID]
      
        allElems = table(cs[[group_by]])
        clusterList = list()
        clusterList[["cluster"]] = clusterID
        clusterList[["all"]] = nrow(cs[[group_by]])
        cs = NULL
        
        for (grp in names(allElems))
        {
          if (!is.null(prefix))
          {
            ngrp = paste(prefix, grp, sep=".")
          } else {
            ngrp = grp
          }
          clusterList[[ngrp]] = allElems[[grp]]
        }
    }
    
    
    
    cellCounts[[clusterID]] = clusterList
    
  }
  df_bysamplerep = cellCounts %>% map(as.data.frame) %>% bind_rows()
  df_bysamplerep[is.na(df_bysamplerep)] <- 0
  
  rownames(df_bysamplerep) = df_bysamplerep$cluster
  df_bysamplerep$cluster = NULL
  
  
  if (relative)
  {
    df_bysamplerep = sweep(df_bysamplerep,2,colSums(df_bysamplerep),"/")
    
    if (show.percent)
    {
      df_bysamplerep = df_bysamplerep*100;
    }
  }
  
  totals=t(colSums(df_bysamplerep))
  totals.df = data.frame(totals)
  rownames(totals.df) = "Total"
  df_bysamplerep=rbind(df_bysamplerep, totals.df)
  
  df_bysamplerep = cbind("cluster"=rownames(df_bysamplerep), df_bysamplerep)
  rownames(df_bysamplerep) = NULL
  
  if (!is.null(outname))
  {
    write.table(df_bysamplerep, file=outname, row.names = F,  quote=FALSE, sep='\t')
  }
  
  return(df_bysamplerep)
  #  
}
countByIACellType = getCellCountDF(obj.integrated, prefix="", select_by = "orig_celltype", group_by="pneu_state", relative=F, show.percent=F, outname="count_by_iacelltype.tsv")











