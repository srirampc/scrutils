library(Seurat)


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

integrate_data = function(athaliana.sobj, ndims = 1:30){
    athaliana.anchors = FindIntegrationAnchors(
        object.list = athaliana.sobj, dims = ndims)
    IntegrateData(anchorset = athaliana.anchors, dims = ndims)
}

load_data = function(base.dir, dir.paths, data.names) {
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

    athaliana.sobjects = lapply(1:length(athaliana.list),
        function(i){
            sobj = CreateSeuratObject(counts = athaliana.list[[i]],
                project = data.names[i],
                min.cells = 3,
                min.features = 200)
            sobj = NormalizeData(sobj)
            sobj = FindVariableFeatures(sobj, selection.method="vst",
                    nfeatures=2000)
            sobj
        })
    cat("Objects loaded and Normalized")
    athaliana.sobjects
}

athaliana.sobj = load_data(root.dir, control.dir.paths,
                            project.names)
athaliana.integrated = integrate_data(athaliana.sobj)

