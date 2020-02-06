library(Seurat)
library(cowplot)
library(ggplot2)
library(harmony)

root.dir = "/data/scRNAseq/athaliana/"
root.dir = "/nv/hswarm1/schockalingam6/data2/scRNAseq/data/"
root.dir = "/project/schockalingam6/scRNAseq/data/athaliana/"
control.dir.paths = c(
"E-CURD-5/control",
"E-GEOD-123013/control",
"E-GEOD-121619/control",
"E-CURD-4/control",
"E-ENAD-30/control"
)

project.names = c(
"E-CURD-5",
"E-GEOD-123013",
"E-GEOD-121619",
"E-CURD-4",
"E-ENAD-30"
)


plot.names = c(
"CRD5",
"GEO3",
"GEO6",
"CRD4",
"NAD3"
)

load_matrices = function(base.dir, dir.paths, data.names) {

    athaliana.list = lapply(1:length(dir.paths),
        function(i){
            mtx = Read10X(paste(base.dir, dir.paths[i], sep=""))
            mtx
    })
    names(athaliana.list) = data.names
    print(sapply(athaliana.list, dim))

    common.gene.names = rownames(athaliana.list[[1]])
    for(i in 2:length(athaliana.list)){
    common.gene.names = intersect(common.gene.names,
                rownames(athaliana.list[[i]]))
    }

    for(i in 1:length(athaliana.list)){
        athaliana.list[[i]] = athaliana.list[[i]][common.gene.names,]
    }
    print(sapply(athaliana.list, dim))
    athaliana.list
}

athaliana.mlist = load_matrices(root.dir, control.dir.paths,
                            project.names)
#ctx =  cbind(stim.sparse, ctrl.sparse)
athaliana.combmat = do.call("cbind", athaliana.mlist)
athaliana = CreateSeuratObject(counts = athaliana.combmat, project = "ATHSC", min.cells = 5)
athaliana = athaliana %>% Seurat::NormalizeData(verbose = FALSE)
athaliana = athaliana %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
athaliana = athaliana %>% ScaleData(verbose = FALSE) 
#
athaliana = athaliana %>% RunPCA(pc.genes = athaliana@var.genes, npcs = 20, verbose = FALSE)
#
#
athaliana.stim = lapply(1:length(athaliana.mlist), 
			function(i){
				c(rep(plot.names[i], ncol(athaliana.mlist[[i]])))
			})
astim = do.call("c", athaliana.stim)
athaliana@meta.data$stim <- astim
#
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = athaliana, reduction = "pca", pt.size = .1, group.by = "stim") #, do.return = TRUE)
p2 <- VlnPlot(object = athaliana, features = "PC_1", group.by = "stim", pt.size = .1) #, do.return = TRUE)
p3 = plot_grid(p1,p2)
ggsave("p3-sath1.pdf", p3, width=12, height=5)
#
#
athaliana <- athaliana %>% RunHarmony("stim")
harmony_embeddings <- Embeddings(athaliana, 'harmony')
#
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = athaliana, reduction = "harmony", pt.size = .1, group.by = "stim") #, do.return = TRUE)
p2 <- VlnPlot(object = athaliana, features = "harmony_1", group.by = "stim", pt.size=.1) #do.return = TRUE, pt.size = .1)
p3 = plot_grid(p1,p2)
#
ggsave("p3-sath2.pdf", p3, width=12, height=5)
#
athaliana <- athaliana %>% RunUMAP(reduction = "harmony", dims = 1:20)
athaliana = athaliana %>% FindNeighbors(reduction = "harmony", dims = 1:20)
athaliana = athaliana %>% FindClusters(resolution = 0.5)
athaliana = athaliana %>% identity()
#
options(repr.plot.height = 4, repr.plot.width = 10)
p3 = DimPlot(athaliana, reduction = "umap", group.by = "stim", pt.size = .1, split.by = 'stim')
ggsave("p3-aths3.pdf", p3, width=10, height=4)
#
options(repr.plot.height = 4, repr.plot.width = 6)
p3 = DimPlot(athaliana, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("p3-aths4.pdf", p3, width=6, height=4)
