library(argparser, quietly=TRUE)
library(Seurat, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(sva, quietly=TRUE)
source("data_utils.R")
source("qc_utils.R")
source("plot_utils.R")

combined_umap = function(athaliana, out.dir, reduce_by="pca",
             dims=1:30, img.option="png"){
    #
    athaliana = cluster_umap_seurat(athaliana, reduce_by, 0.5, dims)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="umap",group="dataset", split="dataset",
         width=10, height=4,
         out_file=paste(out.dir, 
                paste(reduce_by, "-combat-umap-grouped.", 
                  img.option, sep=""), 
                sep="/"))
    #
    dim_plot(athaliana, reduce_by="umap",group=NULL, label=TRUE,
         width=6, height=4,
         out_file=paste(out.dir, 
                paste(reduce_by, "-combat-umap-integrated.", 
                  img.option, sep=""), 
                sep="/"))
    athaliana
}


combined_tsne = function(athaliana, out.dir, reduce_by="pca", 
             dims=1:30, img.option="png"){
    #
    athaliana = cluster_tsne_seurat(athaliana, reduce_by, 0.5, dims)
    #
    print(athaliana@reductions)
    dim_plot(athaliana, reduce_by="tsne",group="dataset", split="dataset",
         width=10, height=4,
         out_file=paste(out.dir, 
                paste(reduce_by, "-combat-tsne-grouped.",
                  img.option, sep=""), 
                sep="/"))
    #
    dim_plot(athaliana, reduce_by="tsne",group=NULL, label=TRUE,
         width=6, height=4,
         out_file=paste(out.dir, 
                paste(reduce_by, "-combat-tsne-integrated.", 
                  img.option, sep=""), 
                sep="/"))
    athaliana
}


combat_cluster = function(root.dir, data.file, out.dir, 
                        qc.flag, vis.option, img.option,
                        inc_list_file, exec_list_file,
                        gen_markers, dot_markers_file){
    data.df = read.csv(data.file, header=TRUE, stringsAsFactors=FALSE)
    print(head(data.df))
    expt.dir.paths = data.df[,'dir.paths']
    short.names = data.df[,'short.names']
    project.names = data.df[, 'project.names']
    
    athaliana.mlist = if(as.logical(qc.flag)){
        qcload_10X_matrices(root.dir, expt.dir.paths,
                            project.names, qc_normalize_matrix,
                            inc_list_file, exec_list_file)
    } else {
        load_10X_matrices(root.dir, expt.dir.paths,
                           project.names)
    }
    clst = combined_expr_matrix(athaliana.mlist, short.names)
    athaliana.combmat = clst[[1]]
    athaliana.batch_names = clst[[2]]
    # Use combat with parametric estimation, no plots
    athaliana.combmat = ComBat(dat=athaliana.combmat, 
                               batch=athaliana.batch_names, 
                               mod=NULL, par.prior=TRUE, prior.plots=FALSE)
    athaliana.integrated = matrix_seurat_object(athaliana.combmat, 
                                                athaliana.batch_names, 
                                                project = "ATHSC", npcs=30)
    if(vis.option == "umap"){
        athaliana.integrated = combined_umap(athaliana.integrated, out.dir,
                         "pca", 1:30, img.option)
    }
    if(vis.option == "tsne"){
        athaliana.integrated = combined_tsne(athaliana.integrated, out.dir, 
                         "pca", 1:30, img.option)
    }
    if(gen_markers) {
        mkdf = FindAllMarkers(athaliana.integrated)
        write.table(mkdf, paste(out.dir, paste(vis.option, "-seurat-markers.tsv", 
                                sep=""),
                    sep="/"), 
                    row.names=FALSE, sep="\t")
    }

    if(!is.na(dot_markers_file)){
        mdf = read.table(dot_markers_file, stringsAsFactors=F, sep="\t")
        for(dxf in mdf$V1){
          gdf = read.table(dxf, header=T, 
                           sep="\t", stringsAsFactors=F)
          genes.plot = gdf$ID
          pfx = gsub("\\.tsv", "" , basename(dxf))
          npresent = genes.plot %in% rownames(athaliana)
          cat(dxf, pfx, sum(npresent), sum(!npresent), "\n")
          if(sum(npresent) > 0) {
             dot_fname = paste(out.dir, 
                            paste(pfx, "dot-markers-seurat.", 
                               img.option, sep=""), 
                         sep="/")
             px = DotPlot(athaliana.integrated, features=genes.plot) +
                      ggtitle(pfx)
             ggsave(dot_fname, px, width=10, height=4)
          }
        }
    }

    athaliana.integrated
}
# Create a parser
p <- arg_parser("Pre-process, Normalize, Integrate w. Combat and Cluster")

# Add command line arguments
p <- add_argument(p, "root_dir", help="Root directory of datasets", type="character")
p <- add_argument(p, "data_file_csv", 
                  help="CSV file with dataset info (See ath.control.csv for example)",
                  type="character")
p <- add_argument(p, "out_dir", help="Output directory", type="character")
p <- add_argument(p, "--qc", help="Flag to indicate to preform qc", short='-q', default=TRUE)
p <- add_argument(p, "--vis", help="Visualization option should be tsne/umap ", short='-v', default='umap')
p <- add_argument(p, "--img", help="Output image option should be one png/pdf ", short='-g', default='png')
p <- add_argument(p, "--exc [glist]", help="File containing list of inc. genes (default:None)", short='-e', default=NULL)
p <- add_argument(p, "--inc [glist]", help="File containing list of exc. genes (default:None)", short='-i', default=NULL)
p <- add_argument(p, "--gen_markers", help="Flag to indicate to genreate all markers", short='-m', default=FALSE)
p <- add_argument(p, "--dot_markers [mfile]", help="Generate marker dot plots for each list of markers list in mfile. Each line in mfile is location for a list of markers", short='-d', default=NULL)

# Parse the command line arguments
argv <- parse_args(p)

if((argv$vis == "tsne" || argv$vis == "umap") && 
   (argv$img == "png" || argv$img == "pdf")) {
         combat_cluster(argv$root_dir, argv$data_file_csv,
                argv$out_dir, argv$qc, argv$vis, argv$img, 
                argv$inc, argv$exec, 
                argv$gen_markers, argv$dot_markers)
    
} else {
    print("Invalid image/visualization option.")
    print(p)
}

