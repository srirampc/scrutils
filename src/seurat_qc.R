library(Seurat)
library(ggplot2)
main = function(data.dir, plot.prefix){
scr.data = Read10X(data.dir = data.dir)
scrj = CreateSeuratObject(counts = scr.data, project = "SCRNA",
			     min.cells = 3, min.features = 200)
scrj[["percent.mt"]] = PercentageFeatureSet(scrj, pattern = "^ATMG")
p1 = VlnPlot(scrj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p3 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p4 = CombinePlots(plots = list(p2, p3))
ggsave(paste(plot.prefix, "-qcvln.pdf", sep=""), p1)
ggsave(paste(plot.prefix, "-fscatter.pdf", sep=""), p4, width=14, height=7)

scrj = NormalizeData(scrj, normalization.method = "LogNormalize", scale.factor = 10000)
scrj = FindVariableFeatures(scrj, selection.method = "vst", nfeatures = 2000)
top10 = head(VariableFeatures(scrj), 10)
p5 = VariableFeaturePlot(scrj)
p6 = LabelPoints(plot = p5, points = top10, repel = TRUE)
p7 = CombinePlots(plots = list(p5, p6))
ggsave(paste(plot.prefix, "-vfeat.pdf", sep=""), p7, width=14, height=7)
all.genes = rownames(scrj)
scrj = ScaleData(scrj) #, features = all.genes)
paste("Data Scaled")
scrj = RunPCA(scrj, features = VariableFeatures(object = scrj))
p8 = DimPlot(scrj, reduction = "pca")
ggsave(paste(plot.prefix, "-pca.pdf", sep=""), p8)
}


args = commandArgs(trailingOnly=TRUE)
if(length(args) == 2){
   main(args[1], args[2])
} else {
   cat("Need 2 args : data.dir, plot.prefix")
}
