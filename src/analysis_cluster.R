library(argparser, quietly=TRUE)
library(tidyr)
library(dplyr)

union_markers_analysis_top2 <- function(known_mg_file, predict_mg_file, out_file) {
    vx = read.table(known_mg_file, header=T, stringsAsFactors=F)
    tx = read.table(predict_mg_file, header=T, stringsAsFactors=F)
    cxnames = unique(tx$cluster)
    vxnames = unique(vx$Type)

    type_fx = Vectorize(
        function(x){ 
            nrow(vx[vx$Type == as.character(x), ])
        })
    cluster_fx = Vectorize(
        function(x){ 
            nrow(tx[tx$cluster == as.numeric(x), ] )
        })
    common_fx = Vectorize(
        function(x, y){
            vxgenes = vx[vx$Type == x, "Gene"] ; 
            tx2 = tx[tx$cluster == y,] 
            tx2 = tx2[order(tx2$p_val_adj), ]; 
            maxl = min(500, nrow(tx2));  
            cxgenes = tx2[1:maxl, "gene"];
            common_genes = intersect(vxgenes, cxgenes)
            length(common_genes)
        })
    common_fx_top2 = Vectorize(
        function(x, y, z){
            vxgenes = vx[vx$Type == x, "Gene"];
            
            ty = tx[tx$cluster == y,] 
            ty = ty[order(ty$p_val_adj), ]; 
            maxly = min(500, nrow(ty));  
            cxgenesy = ty[1:maxly, "gene"];
            
            tz = tx[tx$cluster == z,] 
            tz = tz[order(tz$p_val_adj), ]; 
            maxlz = min(500, nrow(tz));  
            cxgenesz = tz[1:maxlz, "gene"];
            
            cxgenes = union(cxgenesy, cxgenesz)
            length(intersect(vxgenes, cxgenes))
        })

     
    adf = crossing(data.frame(Type=vxnames), 
                   data.frame(Cluster1=cxnames),
                   data.frame(Cluster2=cxnames)) %>% 
            mutate(TypeSize = type_fx(Type)) %>% 
            mutate(Cluster1Size = cluster_fx(Cluster1)) %>%
            mutate(Cluster2Size = cluster_fx(Cluster2)) %>%
            mutate(Cluster1Common = common_fx(Type, Cluster1)) %>%
            mutate(Cluster2Common = common_fx(Type, Cluster2)) %>%
            mutate(CommonGenes = common_fx_top2(Type, Cluster1, Cluster2))
   #head(adf)
   #write.table(as.data.frame(adf),  "test.tsv", sep="\t", )
   bdf = adf %>% group_by(Type) %>%
            summarize(
               KnwonMarkersSize2 = max(TypeSize),
               MaxCommon2 = max(CommonGenes),
               NoMaxClusters2 = length(
                  which(CommonGenes == max(CommonGenes))),
               MaxClusters1 = paste(
                  Cluster1[which(CommonGenes == max(CommonGenes))], 
                  collapse = ", "),
               MaxClusters2 = paste(
                  Cluster2[which(CommonGenes == max(CommonGenes))], 
                  collapse = ", ")
            )
   #head(bdf)
   bdf
}


max_markers_analysis <- function(known_mg_file, predict_mg_file, out_file) {
    vx = read.table(known_mg_file, 
                    # "../data/mg-all-genes-tidy.tsv", 
                    header=T, stringsAsFactors=F)
    tx = read.table(predict_mg_file, 
                    # "athaliana/analysis/control/umap-seurat-markers.tsv", 
                    header=T, stringsAsFactors=F)
    cxnames = unique(tx$cluster)
    vxnames = unique(vx$Type)

    type_fx = Vectorize(
        function(x){ 
            nrow(vx[vx$Type == as.character(x), ] )
        })
    cluster_fx = Vectorize(
        function(x){ 
            nrow(tx[tx$cluster == as.numeric(x), ] )
        })
    common_fx = Vectorize(
        function(x, y){
            vxgenes = vx[vx$Type == x, "Gene"] ; 
            tx2 = tx[tx$cluster == y,] 
            tx2 = tx2[order(tx2$p_val_adj), ]; 
            maxl = min(500, nrow(tx2));  
            cxgenes = tx2[1:maxl, "gene"];
            common_genes = intersect(vxgenes, cxgenes)
            length(common_genes)
        })

#    for(i in 1:length(vxnames)) { 
#        vxgenes = vx[vx$Type == vxnames[i], "Gene"] ; 
#        for(j in 1:length(cxnames)) { 
#            #cxgenes = tx[tx$cluster == cxnames[j], "gene"]; 
#            tx2 = tx[tx$cluster == cxnames[j],] 
#            tx2 = tx2[order( tx2$p_val_adj), ]; 
#            maxl = min(500, nrow(tx2));  
#            cxgenes = tx2[1:maxl, "gene"];
#            common_genes = intersect(vxgenes, cxgenes)
#            cat(vxnames[i], cxnames[j], 
#                length(vxgenes), length(cxgenes), 
#                length(common_genes), "\n")
#        }
#    }
     adf = crossing(data.frame(Type=vxnames), data.frame(Cluster=cxnames)) %>% 
            mutate(TypeSize = type_fx(Type)) %>% 
            mutate(ClusterSize = cluster_fx(Cluster)) %>%
            mutate(CommonGenes = common_fx(Type, Cluster))
    # head(adf)
    bdf = adf %>% group_by(Type) %>%
            summarize(
               KnwonMarkersSize = max(TypeSize),
               MaxCommon = max(CommonGenes),
               MaxClusters = paste(
                  Cluster[which(CommonGenes == max(CommonGenes))], 
                  collapse = ", "),
               ClusterSizes = paste(
                  ClusterSize[which(CommonGenes == max(CommonGenes))], 
                  collapse = ", ")
             )
    bdf 
}
main <- function(known_mg_file, predict_mg_file, out_file) {
   bdf1 = max_markers_analysis(known_mg_file, predict_mg_file, 
                               out_file)
   bdf2 = union_markers_analysis_top2(known_mg_file, 
                                      predict_mg_file, 
                                      out_file)
   bdf = bdf1 %>% inner_join(bdf2)
   head(bdf)
   write.table(as.data.frame(bdf),  out_file, sep="\t")
}

# Create a parser
p <- arg_parser("Percentage Genes")
p <- add_argument(p, "known_mg_file", 
                  help="Tidy format of all the known marker genes", 
                  type="character")
p <- add_argument(p, "predict_mg_file", 
                  help="Predicted marker genes", 
                  type="character")
p <- add_argument(p, "out_file", 
                  help="Output file", 
                  type="character")

# Parse the command line arguments
argv <- parse_args(p)
main(argv$known_mg_file, argv$predict_mg_file, argv$out_file)
