library(Seurat)
library(ggplot2)

seurat_allqc_plot = function(data.dir,plot.prefix, plot.suffix){
    scr.data = Read10X(data.dir = data.dir)
    scrj = CreateSeuratObject(counts = scr.data, project = "SCRNA",
                    min.cells = 3, min.features = 200)
    scrj[["percent.mt"]] = PercentageFeatureSet(scrj, pattern = "^ATMG")
    p1 = VlnPlot(scrj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    p2 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    p3 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p4 = CombinePlots(plots = list(p2, p3))
    ggsave(paste(plot.prefix, "/qcvln", plot.suffix, sep=""), p1)
    ggsave(paste(plot.prefix, "/fscatter", plot.suffix, sep=""), p4, width=14, height=7)

    scrj = NormalizeData(scrj, normalization.method = "LogNormalize", scale.factor = 10000)
    scrj = FindVariableFeatures(scrj, selection.method = "vst", nfeatures = 2000)
    top10 = head(VariableFeatures(scrj), 10)
    p5 = VariableFeaturePlot(scrj)
    p6 = LabelPoints(plot = p5, points = top10, repel = TRUE)
    p7 = CombinePlots(plots = list(p5, p6))
    ggsave(paste(plot.prefix, "/vfeat", plot.suffix, sep=""), p7, width=14, height=7)
    all.genes = rownames(scrj)
    scrj = ScaleData(scrj) #, features = all.genes)
    paste("Data Scaled")
    scrj = RunPCA(scrj, features = VariableFeatures(object = scrj))
    p8 = DimPlot(scrj, reduction = "pca")
    ggsave(paste(plot.prefix, "/pca", plot.suffix, sep=""), p8)
}


seurat_ba_qcplot = function(data.dir, all_thrs, dirx, out.dir, out.prefix, plot.suffix){
    scr.data = Read10X(data.dir = data.dir)
    scrj = CreateSeuratObject(counts = scr.data, project = data.dir)
    scrj[["percent.mcg"]] = PercentageFeatureSet(scrj, pattern = "^AT[MC]G")

    p2 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "percent.mcg")
    p3 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p4 = CombinePlots(plots = list(p2, p3))
    fname = paste(out.dir, dirx, paste(out.prefix, "-seurat-before-qc.png", sep=""), 
		  sep="/"  )
    ggsave(fname, p4, width=14, height=7)

    mcg_threshold = all_thrs[1]
    min_cells = all_thrs[2]
    min_features = all_thrs[3]
    cat("MCG minc", data.dir, mcg_threshold, min_cells, "\n")
    scrj = CreateSeuratObject(counts = scr.data, project = data.dir,
			      min.cells = all_thrs[2], min.features = min_features)
    scrj[["percent.mcg"]] = PercentageFeatureSet(scrj, pattern = "^AT[MC]G")
    #srch =  subset(pbmc, subset = nCount_RNA > 200 & nFeature_RNA < 2500 & percent.mcg < mcg_threshold)
    #scrj =  subset(scrj, subset = nCount_RNA > 200 & percent.mcg < mcg_threshold)
    expr1 = FetchData(object=scrj, vars="nFeature_RNA")
    scrj1 = scrj[, expr1 > 200]
    expr2 = FetchData(object=scrj1, vars="percent.mcg")
    scrj2 = scrj1[, expr2 < mcg_threshold]
    print(scrj2)

    p2 = FeatureScatter(scrj2, feature1 = "nCount_RNA", feature2 = "percent.mcg")
    p3 = FeatureScatter(scrj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p4 = CombinePlots(plots = list(p2, p3))
    fname = paste(out.dir, dirx, paste(out.prefix, "-seuerat-after-qc.png", sep=""), 
		  sep="/"  )
    ggsave(fname, p4, width=14, height=7)
}

