library(Seurat)
source("data_utils.R")

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



athaliana.sobj = load_10X_seurat_objects(root.dir, control.dir.paths,
                            project.names)
athaliana.integrated = integrate_seurat_objects(athaliana.sobj)

