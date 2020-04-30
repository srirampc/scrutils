library(Seurat)
source("data_utils.R")


seurat_cluster = function(root.dir, data.file.csv){
    data.df = read.csv(data.file, header=TRUE)
    expt.dir.paths = data.df$dir.paths
    short.names = data.df$short.names
    project.names = data.df$project.names

    athaliana.sobj = load_10X_seurat_objects(root.dir, expt.dir.paths,
                                project.names)
    athaliana.integrated = integrate_seurat_objects(athaliana.sobj)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 2){
    seurat_cluster(args[1], args[2])
}  else {
    print(args)
    print("Usage: Rscript seurat_cluster.R root_dir data.file.csv")
}
