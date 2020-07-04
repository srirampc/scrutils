
library(Seurat, quietly=TRUE)
library(harmony, quietly=TRUE)
library(stringr, quietly=TRUE)

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

load_10X_matrices_union = function(base.dir, dir.paths, data.names) {
    mat10x.list = lapply(1:length(dir.paths),
        function(i){
            mtx = Read10X(paste(base.dir, dir.paths[i], sep=""))
            mtx
    })
    names(mat10x.list) = data.names
    print(sapply(mat10x.list, dim))

    all.gene.names = rownames(mat10x.list[[1]])
    for(i in 2:length(mat10x.list)){
         all.gene.names = union(all.gene.names,
                                rownames(mat10x.list[[i]]))
    }

    for(i in 1:length(mat10x.list)){
        missing.genes = setdiff(all.gene.names, 
                                rownames(mat10x.list[[i]]))
        mtx1 = as.matrix(mat10x.list[[i]])
	mtx2 = matrix(rep(0, dim(mtx1[2]) * length(missing.genes)), 
                      nrow = length(missing.genes),
                      ncol = dim(mtx1[2]))
	rownames(mtx2) = missing.genes
	mtx = rbind(mtx1, mtx2)
        mat10x.list[[i]] = mtx[all.gene.names,]
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


integrate_seurat_objects = function(athaliana.sobj, anc.feat, ndims = 1:30){
    athaliana.anchors = FindIntegrationAnchors(
        object.list = athaliana.sobj, anchor.features=anc.feat,
        dims = ndims)
    IntegrateData(anchorset = athaliana.anchors, dims = ndims)
}


combined_expr_matrix = function(data.mlist, short.names){
    data.batch_names = lapply(1:length(data.mlist), 
                function(i){
                    c(rep(short.names[i], ncol(data.mlist[[i]])))
                })
    data.unames = lapply(1:length(data.mlist), 
                function(i){
                    str_c(short.names[i], 1:ncol(data.mlist[[i]]), sep="-")
                })
    data.batch_names = do.call("c", data.batch_names)
    data.unames = do.call("c", data.unames)
    data.combmat = do.call("cbind", data.mlist)
    rownames(data.combmat) = rownames(data.mlist[[1]])
    colnames(data.combmat) = data.unames
    list(data.combmat, data.batch_names)
}

matrix_seurat_object = function(data.combmat, data.batch_names,
                                normalize = TRUE,
                                project = "ATHSC", dim_name="dataset", 
                                mincells=5, nfeat=2000, npcs=20){
    data.sobj = CreateSeuratObject(counts = data.combmat, project = project, 
		                   min.cells = mincells)
    data.sobj@meta.data[dim_name] <- data.batch_names
    if(normalize == TRUE){
        data.sobj = data.sobj %>% Seurat::NormalizeData(verbose = FALSE)
    }
    data.sobj = data.sobj %>% ScaleData(verbose = FALSE)
    data.sobj = data.sobj %>% FindVariableFeatures(selection.method = "vst", 
                                                   nfeatures = nfeat,
                                                   verbose=FALSE)
    #
    if(npcs > 0) {
        data.sobj = data.sobj %>% RunPCA(pc.genes = (data.sobj)@var.genes, 
                                         npcs = npcs, verbose = FALSE)
    }
    data.sobj
}

combined_seurat_object = function(data.mlist, short.names, normalize=TRUE,
                                  project = "ATHSC", dim_name="dataset",
                                  mincells=5, nfeat=2000, npcs=20) {
    clst = combined_expr_matrix(data.mlist, short.names)
    data.combmat = clst[[1]]
    data.batch_names = clst[[2]]
    matrix_seurat_object(data.combmat, data.batch_names, normalize,
                         project, dim_name,
                         mincells, nfeat, npcs)
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

qcload_10X_matrices = function(base.dir, dir.paths, data.names, qc.function,
                               inc_list_file, exc_list_file) {
    mat10x.list = lapply(1:length(dir.paths),
        function(i){
            # dfx = qc_normalize(paste(base.dir, dir.paths[i], sep="/"))
            # ctmtx = counts(dfx)
            # print(dim(ctmtx))
            # rownames(ctmtx) = rownames(dfx)
            # colnames(ctmtx) = 1:dim(dfx)[2]
            # ctmtx
            cat("Loading ", dir.paths[i], "...\n")
            qc.function(paste(base.dir, dir.paths[i], sep="/"),
                        inc_list_file, exc_list_file)
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
    cat("Matrix Objects loaded for all datasets",
        sapply(mat10x.list, dim) ,"\n")

    mat10x.list
}


qcload_10X_matrices_union = function(base.dir, dir.paths, data.names, qc.function,
                                     inc_list_file, exc_list_file) {
    mat10x.list = lapply(1:length(dir.paths),
        function(i){
            # dfx = qc_normalize(paste(base.dir, dir.paths[i], sep="/"))
            # ctmtx = counts(dfx)
            # print(dim(ctmtx))
            # rownames(ctmtx) = rownames(dfx)
            # colnames(ctmtx) = 1:dim(dfx)[2]
            # ctmtx
            cat("Loading ", dir.paths[i], "...\n")
            qc.function(paste(base.dir, dir.paths[i], sep="/"),
                        inc_list_file, exc_list_file)
    })
    names(mat10x.list) = data.names
    print(sapply(mat10x.list, dim))

    all.gene.names = rownames(mat10x.list[[1]])
    for(i in 2:length(mat10x.list)){
        all.gene.names = union(all.gene.names,
                               rownames(mat10x.list[[i]]))
    }

    for(i in 1:length(mat10x.list)){
        missing.genes = setdiff(all.gene.names, 
                                rownames(mat10x.list[[i]]))
        mtx1 = as.matrix(mat10x.list[[i]])
	mtx2 = matrix(rep(0, dim(mtx1)[2] * length(missing.genes)), 
                      nrow = length(missing.genes),
                      ncol = dim(mtx1)[2])
	rownames(mtx2) = missing.genes
	mtx = rbind(mtx1, mtx2)
        mat10x.list[[i]] = mtx[all.gene.names,]
    }
    cat("Matrix Objects loaded for all datasets",
        sapply(mat10x.list, dim) ,"\n")

    mat10x.list
}

qcload_10X_seurat_objects = function(base.dir, dir.paths, data.names, 
                                     short.names, qc.function,
                                     inc_list_file, exc_list_file,
                                     dim.name="dataset", nfeats=2000,
		                     min.cells = 2, min.feats=200) {

    mat10X.sobjects = lapply(1:length(dir.paths),
        function(i){
            cat("Loading ", dir.paths[i], "...\n")
            cmatx = qc.function(paste(base.dir, dir.paths[i], sep="/"),
                                inc_list_file, exc_list_file)
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
    cat("Seurat Objects loaded for all datasets\n")
    mat10X.sobjects
}
