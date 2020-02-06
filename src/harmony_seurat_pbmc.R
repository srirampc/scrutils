library(Seurat)
library(cowplot)
library(harmony)

pbmc_rdata = 'hdata/pbmc_stim.RData'
plot_outdir = "."
load(pbmc_rdata)
ctx =  cbind(stim.sparse, ctrl.sparse)
pbmc = CreateSeuratObject(counts = ctx, project = "PBMC", min.cells = 5)
pbmc = pbmc %>% Seurat::NormalizeData(verbose = FALSE)
pbmc = pbmc %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
pbmc = pbmc %>% ScaleData(verbose = FALSE) 

pbmc = pbmc %>% RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)


pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse)))

options(repr.plot.height = 5, repr.plot.width = 12)
p1 = DimPlot(object = pbmc, reduction = "pca", pt.size = .1, group.by = "stim") #, do.return = TRUE)
p2 = VlnPlot(object = pbmc, features = "PC_1", group.by = "stim", pt.size = .1) #, do.return = TRUE)
p3 = plot_grid(p1,p2)
ggsave(paste(plot_outdir, "p3-s1.pdf", sep=""), p3, height=5, width=12)


pbmc <- pbmc %>% RunHarmony("stim")
harmony_embeddings <- Embeddings(pbmc, 'harmony')

options(repr.plot.height = 5, repr.plot.width = 12)
p1 = DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "stim") #, do.return = TRUE)
p2 = VlnPlot(object = pbmc, features = "harmony_1", group.by = "stim", pt.size=.1) #do.return = TRUE, pt.size = .1)
p3 = plot_grid(p1,p2)
ggsave(paste(plot_outdir, "p3-s2.pdf", sep=""), p3, height=5, width=12)

pbmc = pbmc %>%
    RunUMAP(reduction = "harmony", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    FindClusters(resolution = 0.5) %>%
    identity()

options(repr.plot.height = 4, repr.plot.width = 10)
p3 = DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1, split.by = 'stim')
ggsave(paste(plot_outdir, "p3-s3.pdf", sep=""), p3)

options(repr.plot.height = 4, repr.plot.width = 6)
p3 = DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = .1)
ggsave(paste(plot_outdir, "p3-s4.pdf", sep=""), p3)
