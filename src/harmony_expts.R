library(dplyr)
library(ggplot2)
library(harmony)
library(cowplot)
library(Seurat)

do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE,
                       do_labels = TRUE, nice_names,
                       palette_use = c(`jurkat` = '#810F7C', `t293` = '#D09E2D',`half` = '#006D2C'),
                       pt_size = 4, point_size = .5, base_size = 12,
                       do_points = TRUE, do_density = FALSE, h = 6, w = 8) {
    umap_use <- umap_use[, 1:2]
    colnames(umap_use) <- c('X1', 'X2')
    plt_df <- umap_use %>% data.frame() %>%
        cbind(meta_data) %>%
        dplyr::sample_frac(1L)
    plt_df$given_name <- plt_df[[label_name]]

    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>%
            subset(nice_name != "" & !is.na(nice_name))

        plt_df[[label_name]] <- plt_df$nice_name
    }

    plt <- plt_df %>%
        ggplot(aes_string("X1", "X2", col = label_name, fill = label_name)) +
        theme_test(base_size = base_size) +
        theme(panel.background = element_rect(fill = NA, color = "black")) +
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                        shape = 16, size = 4)),
               alpha = FALSE) +
        scale_color_manual(values = palette_use) +
        scale_fill_manual(values = palette_use) +
        theme(plot.title = element_text(hjust = .5)) +
        labs(x = "PC 1", y = "PC 2")

    if (do_points)
        plt <- plt + geom_point(shape = '.')
    if (do_density)
        plt <- plt + geom_density_2d()


    if (no_guides)
        plt <- plt + guides(col = FALSE, fill = FALSE, alpha = FALSE)

    if (do_labels) {
        data_labels <- plt_df %>%
            dplyr::group_by_(label_name) %>%
            dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>%
            dplyr::ungroup()

        plt <- plt + geom_label(data = data_labels, label.size = NA,
                        aes_string(label = label_name),
                        color = "white", size = pt_size, alpha = 1,
                        segment.size = 0) +
                guides(col = FALSE, fill = FALSE)
    }

    return(plt)
}

harmony_cell_lines = function(){
    out_dir = "./"
    data(cell_lines)
    V <- cell_lines$scaled_pcs
    meta_data <- cell_lines$meta_data
    p1 = do_scatter(V, meta_data, 'dataset') + 
        labs(title = 'Colored by dataset')
    p2 = do_scatter(V, meta_data, 'cell_type') +
        labs(title = 'Colored by cell type')
    p3 = cowplot::plot_grid(p1, p2)
    ggsave(paste(out_dir, "p3-1.pdf", sep=""), p3)
    harmony_embeddings <- harmony::HarmonyMatrix(
        V, meta_data, 'dataset', do_pca = FALSE, verbose=FALSE
    )
    p1 <- do_scatter(harmony_embeddings, meta_data, 'dataset') + 
        labs(title = 'Colored by dataset')

    p2 <- do_scatter(harmony_embeddings, meta_data, 'cell_type') +
        labs(title = 'Colored by cell type')
    p3 = cowplot::plot_grid(p1, p2, nrow = 1)

    ggsave(paste(out_dir, "p3-2.pdf", sep=""), p3)
}

harmony_pbmc = function(){
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

}

harmony_cell_lines()
harmony_pbmc()