import numpy as np
import pandas as pd
import scanpy as sc
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sb
from anndata import AnnData

def scanpy_qc(adata, out_dir, out_prefix, mito_prefix):
    adata.var_names_make_unique()
    #
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith(mito_prefix) for gene in adata.var_names]
    adata.obs['mt_frac'] = adata[:, adata.var['mt']].X.sum(1).A.squeeze()/adata.obs['total_counts']
    fig, axs = plt.subplots(1,4, figsize=(15,4))
    fig.suptitle('Covariates for filtering')
    sb.distplot(adata.obs['total_counts'], kde=False, ax = axs[0])
    sb.distplot(adata.obs['total_counts'][adata.obs['total_counts']<10000], kde=False, bins=40, ax = axs[1])
    sb.distplot(adata.obs['n_genes_by_counts'], kde=False, bins=60, ax = axs[2])
    sb.distplot(adata.obs['n_genes_by_counts'][adata.obs['n_genes_by_counts']<4000], kde=False, bins=60, ax = axs[3])
    fig.savefig(out_dir + "/visum_qc_histo" + out_prefix )
    with open(out_dir + "/visum_qc_thres.txt", 'w') as ofx:
        sc.pp.filter_cells(adata, min_counts = 5000)
        ofx.write(f'Number of cells after min count filter (>5K) : {adata.n_obs}')
        ofx.write("\n")
        sc.pp.filter_cells(adata, max_counts = 35000)
        ofx.write(f'Number of cells after max count filter (<35K) : {adata.n_obs}')
        ofx.write("\n")
        adata = adata[adata.obs['mt_frac'] < 0.2]
        ofx.write(f'Number of cells after MT filter (< 0.2) : {adata.n_obs}')
        ofx.write("\n")
        sc.pp.filter_cells(adata, min_genes = 3000)
        ofx.write(f'Number of cells after gene filter (min. 3K) : {adata.n_obs}')
        ofx.write("\n")
        sc.pp.filter_genes(adata, min_cells=10)
        ofx.write(f'Number of genes after cell filter (min. 10) : {adata.n_vars}')
        ofx.write("\n")

def seurat_wf_plots(adata, out_dir, out_prefix, mito_prefix):
    sc.pl.highest_expr_genes(adata, n_top=20, save=out_prefix, show=False)
    # Filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    #
    mito_genes = adata.var_names.str.startswith(mito_prefix)
    out_prefix = '_seuratwf' + out_prefix
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                 jitter=0.4, multi_panel=True, save=out_prefix, show=False)
    sc.pl.scatter(adata, x='n_counts', y='percent_mito',
                  save="_pmito"+out_prefix, show=False)
    sc.pl.scatter(adata, x='n_counts', y='n_genes',
                 save="_genes_"+out_prefix, show=False)
    # subset & Normalize
    adata = adata[adata.obs.n_genes < 2500, :]
    adata = adata[adata.obs.percent_mito < 0.05, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Freeze state
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, save=out_prefix, show=False)
    #
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['n_counts', 'percent_mito'], copy=True)
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata, save=out_prefix, show=False)
    sc.pl.pca_variance_ratio(adata, log=True, save=out_prefix, show=False)
    # UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    #sc.tl.paga(adata)
    #sc.pl.paga(adata, plot=False)
    sc.tl.umap(adata)
    sc.pl.umap(adata, save=out_prefix, show=False)
    # t-SNE
    sc.tl.tsne(adata)
    sc.pl.tsne(adata, save=out_prefix, show=False)


def recipe_seurat(adata: AnnData, out_dir, out_prefix,
                  log: bool = True, copy: bool = False, ):
    """
    Taken from 
    https://github.com/theislab/scanpy/blob/master/scanpy/preprocessing/_recipes.py
    """
    if copy: adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5, log=not log)
    sc.pl.filter_genes_dispersion(filter_result, log=not log,
                                  save="_seurat" + out_prefix, show=False)
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    if log: sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    return adata if copy else None



def recipe_zheng17(adata: AnnData, out_dir, out_prefix,
                   n_top_genes: int = 2000, log: bool = True,
                   copy: bool = False,):
    """
    Taken from 
    https://github.com/theislab/scanpy/blob/master/scanpy/preprocessing/_recipes.py
    """
    start = sc.logging.info('running recipe zheng17')
    if copy: adata = adata.copy()
    # only consider genes with more than 1 count
    sc.pp.filter_genes(adata, min_counts=1)
    # normalize with total UMI count per cell
    sc.pp.normalize_total(adata, key_added='n_counts_all')
    filter_result = sc.pp.filter_genes_dispersion(
        adata.X, flavor='cell_ranger', n_top_genes=n_top_genes, log=False
    )
    sc.pl.filter_genes_dispersion(filter_result, log=True,
                                  save="_zheng"+ out_prefix,
                                  show=False)
    # actually filter the genes, the following is the inplace version of
    #     adata = adata[:, filter_result.gene_subset]
    adata._inplace_subset_var(filter_result.gene_subset)  # filter genes
    sc.pp.normalize_total(adata)  # renormalize after filtering
    if log: sc.pp.log1p(adata)  # log transform: X = log(X + 1)
    sc.pp.scale(adata)
    sc.logging.info('    finished', time=start)
    return adata if copy else None

def main(src_dir, out_dir, out_prefix, mito_prefix):
    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)
    sc.settings.figdir = out_dir + "/"
    adata = sc.read_10x_mtx(src_dir, var_names='gene_symbols',
                            cache = True)
    adata.var_names_make_unique()
    adata2 = adata.copy()
    adata3 = adata.copy()
    adata4 = adata.copy()
    seurat_wf_plots(adata, out_dir, out_prefix, mito_prefix)
    recipe_seurat(adata2, out_dir, out_prefix)
    recipe_zheng17(adata3, out_dir, out_prefix)
    scanpy_qc(adata4, out_dir, out_prefix, mito_prefix)



if __name__ == "__main__":
    PROG_DESC = """
    Preprocess plots usin scanpy
    """
    PARSER = argparse.ArgumentParser(description=PROG_DESC)
    PARSER.add_argument("src_dir",
                        help="""Input directrory containing 
                                matrix.mtx, genes.tsv and barcodes.tsv files""")
    PARSER.add_argument("output_dir",
                        help="""Output directory to create all the plots""")
    PARSER.add_argument("out_prefix",
                        help="""Prefix for output plot files""")
    PARSER.add_argument("mito_prefix",
                        help="""Prefix for mitochondria genes""")
    ARGS = PARSER.parse_args()
    print(ARGS.src_dir, ARGS.output_dir,
          ARGS.out_prefix, ARGS.mito_prefix)
    main(ARGS.src_dir, ARGS.output_dir,
         ARGS.out_prefix, ARGS.mito_prefix)
