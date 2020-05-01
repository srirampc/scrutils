
library(Seurat)
library(harmony)
library(stringr)

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

load_10X_seurat_objects = function(base.dir, dir.paths, data.names,
                                    short.names, dim.name="dataset", min.cells = 3,
                                    min.feats=200, nfeats=2000) {
    mat10x.list = load_10X_matrices(base.dir, dir.paths, data.names)

    mat10X.sobjects = lapply(1:length(mat10x.list),
        function(i){
            sobj = CreateSeuratObject(counts = mat10x.list[[i]],
                project = data.names[i],
                min.cells = min.cells,
                min.features = min.feats)
            sobj = NormalizeData(sobj)
            sobj = FindVariableFeatures(sobj, selection.method="vst",
                    nfeatures=nfeats)
            sobj@meta.data[dim.name] <- c(rep(short.names[i], ncol(cmatx)))
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


combined_seurat_object = function(data.mlist, short.names, dim_name="dataset",
				  nfeat=2000, mincells=5, npcs=20) {
    data.cnames = lapply(1:length(data.mlist), 
                function(i){
                    c(rep(short.names[i], ncol(data.mlist[[i]])))
                })
    data.unames = lapply(1:length(data.mlist), 
                function(i){
                    str_c(short.names[i], 1:ncol(data.mlist[[i]]), sep="-")
                })
    data.cnames = do.call("c", data.cnames)
    data.unames = do.call("c", data.unames)
    data.combmat = do.call("cbind", data.mlist)
    rownames(data.combmat) = rownames(data.mlist[[1]])
    colnames(data.combmat) = data.unames
    data.sobj = CreateSeuratObject(counts = data.combmat, project = "ATHSC", 
				   min.cells = mincells)
    data.sobj = data.sobj %>% Seurat::NormalizeData(verbose = FALSE)
    data.sobj = data.sobj %>% FindVariableFeatures(selection.method = "vst", 
						   nfeatures = nfeat)
    data.sobj = data.sobj %>% ScaleData(verbose = FALSE) 
    #
    data.sobj = data.sobj %>% RunPCA(pc.genes = (data.sobj)@var.genes, 
				     npcs = npcs, verbose = FALSE)
    data.sobj@meta.data[dim_name] <- data.cnames
    data.sobj
}

integrate_data_harmony = function(combo.sobj, dim_name="dataset"){
    combo.sobj <- combo.sobj %>% RunHarmony(dim_name)
    harmony_embeddings <- Embeddings(combo.sobj, 'harmony')
    list(combo.sobj, harmony_embeddings)
} 

cluster_umap_seurat = function(data.obj, reduce_by, resolution=0.5, dims=1:20){
    data.obj = data.obj %>% RunUMAP(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindNeighbors(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindClusters(resolution = resolution)
    data.obj = data.obj %>% identity()
    data.obj
}


scran_normalize = function(dfx){
    clusters <- scran::quickCluster(dfx)
    dfx <- scran::computeSumFactors(dfx, clusters=clusters)
    #print(summary(scater::sizeFactors(dfx)))
    dfx <- scater::logNormCounts(dfx)
    dfx
}

cluster_tsne_seurat = function(data.obj, reduce_by, resolution=0.5, dims=1:20){
    data.obj = data.obj %>% RunTSNE(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindNeighbors(reduction = reduce_by, dims = dims)
    data.obj = data.obj %>% FindClusters(resolution = resolution)
    data.obj = data.obj %>% identity()
    data.obj
}

qcload_10X_matrices = function(base.dir, dir.paths, data.names, qc.function) {
    mat10x.list = lapply(1:length(dir.paths),
        function(i){
            # dfx = qc_normalize(paste(base.dir, dir.paths[i], sep="/"))
            # ctmtx = counts(dfx)
            # print(dim(ctmtx))
            # rownames(ctmtx) = rownames(dfx)
            # colnames(ctmtx) = 1:dim(dfx)[2]
            # ctmtx
            qc.function(paste(base.dir, dir.paths[i], sep="/"))
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


qcload_10X_seurat_objects = function(base.dir, dir.paths, data.names, 
                                     short.names, qc.function, dim.name="dataset",
				     min.cells = 3, min.feats=200, nfeats=2000) {

    mat10X.sobjects = lapply(1:length(dir.paths),
        function(i){
            cmatx = qc.function(paste(base.dir, dir.paths[i], sep="/"))
            colnames(cmatx) = str_c(short.names[i], 1:dim(cmatx)[2], sep="_")
            sobj = CreateSeuratObject(counts = cmatx,
                project = data.names[i],
                min.cells = min.cells,
                min.features = min.feats)
            # sobj = NormalizeData(sobj)
            sobj = FindVariableFeatures(sobj, selection.method="vst",
                    nfeatures=nfeats)
            sobj@meta.data[dim.name] <- c(rep(short.names[i], ncol(cmatx)))
            sobj
        })
    cat("Objects loaded")
    mat10X.sobjects
}
