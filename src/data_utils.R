
library("Seurat")
load_10X_matrices = function(base.dir, dir.paths, data.names) {
    mat10x.list = lapply(1:length(dir.paths),
        function(i){
            mtx = Read10X(paste(base.dir, dir.paths[i], sep=""))
            mtx
    })
    names(mat10x.list) = data.names
    print(sapply(mat10x.list, dim))

    common.gene.names = rownames(mat10x.list[[1]])
    for(i in 2:length(mat10x.list)){
    common.gene.names = intersect(common.gene.names,
                rownames(mat10x.list[[i]]))
    }

    for(i in 1:length(mat10x.list)){
        mat10x.list[[i]] = mat10x.list[[i]][common.gene.names,]
    }
    print(sapply(mat10x.list, dim))
    mat10x.list
}

load_10X_seurat_objects = function(base.dir, dir.paths, data.names) {
    mat10x.list = load_10X_matrices(base.dir, dir.paths, data.names)

    mat10X.sobjects = lapply(1:length(mat10x.list),
        function(i){
            sobj = CreateSeuratObject(counts = mat10x.list[[i]],
                project = data.names[i],
                min.cells = 3,
                min.features = 200)
            sobj = NormalizeData(sobj)
            sobj = FindVariableFeatures(sobj, selection.method="vst",
                    nfeatures=2000)
            sobj
        })
    cat("Objects loaded and Normalized")
    mat10X.sobjects
}



integrate_seurat_objects = function(athaliana.sobj, ndims = 1:30){
    athaliana.anchors = FindIntegrationAnchors(
        object.list = athaliana.sobj, dims = ndims)
    IntegrateData(anchorset = athaliana.anchors, dims = ndims)
}


combined_seurat_object = function(data.mlist, short.names) {
    data.combmat = do.call("cbind", data.mlist)
    data.sobj = CreateSeuratObject(counts = data.combmat, project = "ATHSC", min.cells = 5)
    data.sobj = data.sobj %>% Seurat::NormalizeData(verbose = FALSE)
    data.sobj = data.sobj %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    data.sobj = data.sobj %>% ScaleData(verbose = FALSE) 
    #
    data.sobj = data.sobj %>% RunPCA(pc.genes = (data.sobj)@var.genes, npcs = 20, verbose = FALSE)
    data.cnames = lapply(1:length(data.mlist), 
                function(i){
                    c(rep(short.names[i], ncol(data.mlist[[i]])))
                })
    data.cnames = do.call("c", data.cnames)
    data.sobj@meta.data$stim <- data.cnames
    data.sobj
}

integrate_data_harmony = function(combo.sobj){
    library(harmony)
    combo.sobj <- combo.sobj %>% RunHarmony("stim")
    harmony_embeddings <- Embeddings(combo.sobj, 'harmony')
    list(comb.sobj, harmony_embeddings)
} 

cluster_umap_seurat = function(data.obj, reduce_by, resolution=0.5, dims=1:20){
    library(Seurat)
    data.obj = data.obj %>% RunUMAP(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindNeighbors(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindClusters(resolution = resolution)
    data.obj = data.obj %>% identity()
    data.obj
}


scran_normalize = function(dfx){
    clusters <- scran::quickCluster(dfx)
    dfx <- scran::computeSumFactors(dfx, clusters=clusters)
    summary(scran::sizeFactors(dfx))
    dfx <- scater::logNormCounts(dfx)
    dfx
}

