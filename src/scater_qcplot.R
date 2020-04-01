library(scater)
library("DropletUtils")
library(Seurat)
library(ggplot2)
seurat_qcplot = function(data.dir, all_thrs, out.dir, out.prefix, plot.suffix){
    scr.data = Read10X(data.dir = data.dir)
    scrj = CreateSeuratObject(counts = scr.data, project = data.dir)
    scrj[["percent.mcg"]] = PercentageFeatureSet(scrj, pattern = "^AT[MC]G")

    p2 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "percent.mcg")
    p3 = FeatureScatter(scrj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p4 = CombinePlots(plots = list(p2, p3))
    fname = paste(out.dir, data.dir, paste(out.prefix, "-seurat-before-qc.png", sep=""), 
		  sep="/"  )
    ggsave(fname, p4, width=14, height=7)

    mcg_threshold = all_thrs[1]
    min_cells = all_thrs[2]
    cat("MCG minc", data.dir, mcg_threshold, min_cells, "\n")
    scrj = CreateSeuratObject(counts = scr.data, project = data.dir,
			      min.cells = all_thrs, min.features = 200)
    scrj[["percent.mcg"]] = PercentageFeatureSet(scrj, pattern = "^AT[MC]G")
    #srch =  subset(pbmc, subset = nCount_RNA > 200 & nFeature_RNA < 2500 & percent.mcg < mcg_threshold)
    #scrj =  subset(scrj, subset = nCount_RNA > 200 & percent.mcg < mcg_threshold)
    expr1 = FetchData(object=scrj, vars="nFeature_RNA")
    scrj1 = scrj[, expr1 > 200]
    expr2 = FetchData(object=scrj1, vars="percent.mcg")
    scrj2 = scrj1[, expr2 < mcg_threshold]
    print(scrj2)

    p2 = FeatureScatter(scrj2, feature1 = "nCount_RNA", feature2 = "percent.mcg")
    p3 = FeatureScatter(scrj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p4 = CombinePlots(plots = list(p2, p3))
    fname = paste(out.dir, data.dir, paste(out.prefix, "-seuerat-after-qc.png", sep=""), 
		  sep="/"  )
    ggsave(fname, p4, width=14, height=7)
}

plot_feat_counts_flip_hist = function(per.feat, dfx, dirx, fname){
    featcts = per.feat$detected * dim(dfx)[2]/100.0
    log_featcts = log10(featcts)
    png(file=fname)
    xhist = hist(log_featcts, breaks=200,plot = FALSE)
    plot(x = xhist$mids, y = xhist$counts, type='h',
            ylab=expression(Log[10]~"No. cells"), xlab="No. genes", main=dirx)
    xrange = (xhist$mids < 2.0) & (xhist$mids > 1.0)  & (xhist$counts > 0)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    #print(xrnmidv)
    abline(v = xrnmidv, col='blue')
    dev.off()
}

plot_counts_flip_hist = function(per.cell, dirx, fname){
    libsize = per.cell$sum
    log_libsize = log10(libsize); 
    png(file=fname)
    xhist = hist(log_libsize, breaks=200,plot = FALSE)
    plot(x = xhist$mids, y = xhist$counts, type='h', xlab="No. cells", ylab=expression(Log[10]~"No. Reads"))
    xrange = (xhist$mids < 2.0) & (xhist$mids > 1.0)  & (xhist$counts > 0)
    #print(sum(xrange))
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    #print(xrnmidv)
    abline(v = xrnmidv, col='blue')
    dev.off()
}

plot_cell_ngenes_flip_hist = function(per.cell, dirx, fname){
    ngenes = per.cell$detected
    png(file=fname)
    xhist = hist(ngenes, breaks=200,plot = FALSE)
    plot(x = xhist$mids, y = xhist$counts, type='h', xlab="No. cells", ylab="No. genes")
    dev.off()

}

plot_mcg_flip_hist = function(per.cell, dirx, fname){
    mcgpct = per.cell$subsets_MCG_percent
    log_mcg = log10(mcgpct)
    png(file=fname)
    xhist = hist(log_mcg, breaks=200,plot = FALSE)
    plot(x = xhist$mids, y = xhist$counts, type='h', xlab="No. cells", ylab=expression(Log[10]~"MCG %"))
    xrange = (xhist$mids < 2.0) & (xhist$mids > 1.0)  & (xhist$counts > 0)
    #print(sum(xrange))
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    #print(xrnmidv)
    abline(v = xrnmidv, col='blue')
    dev.off()

}
plot_feat_avg_flip_hist = function(per.feat, dirx, fname){
    featavg = per.feat$mean 
    log_favg = log10(featavg)
    png(file=fname)
    xhist = hist(log10(featavg), breaks=100,plot = FALSE)
    plot(x = xhist$mids, y = xhist$counts, type='h',
        ylab=expression(Log[10]~"Avg. No. Reads"), xlab="No. cells", 
        main=dirx)
    xrange = (xhist$mids < -1.5) & (xhist$mids > -2.75)  & (xhist$counts > 0)
    #print(sum(xrange))
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    print(xrnmidv)
    abline(v = xrnmidv, col='blue')
    dev.off()
}

plot_flip_hist = function(per.cell, per.feat, dfx, dirx, out_dir, out_prefix){
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-feat-filp-hist.png", sep=""), sep="/"  )
    plot_feat_counts_flip_hist(per.feat, dfx, dirx, fname)
    
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-counts-flip-hist.png", sep=""), sep="/"  )
    plot_counts_flip_hist(per.cell, dirx, fname)
    
    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-ngenes-flip-hist.png", sep=""), sep="/"  )
    plot_cell_ngenes_flip_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-mcgpct-flip-hist.png", sep=""), sep="/"  )
    plot_mcg_flip_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, paste(out_prefix, "-scater-qc-favg-flip-hist.png", sep=""), sep="/"  )
    plot_feat_avg_flip_hist(per.feat, dirx, fname)
}

plot_mcg_hist = function(per.cell, dirx, fname){
    mcgpct = per.cell$subsets_MCG_percent
    log_mcg = log10(mcgpct)
    png(file=fname)
    log_mcgf = log_mcg
    log_mcg = log_mcg[log_mcg > -2.3]
    mcg_nmads = isOutlier(log_mcg, nmads=3, typ="higher")
    mcg_mad = mad(log_mcg)
    mcg_med = median(log_mcg)
    mcg_3mads = mcg_med - 3 * mcg_mad
    mcg_2mads = mcg_med - 2 * mcg_mad
    mcg_1mads = mcg_med - mcg_mad
    mcg_3madsp = mcg_med + 3 * mcg_mad
    mcg_2madsp = mcg_med + 2 * mcg_mad
    mcg_1madsp = mcg_med + mcg_mad
    # cat("3MAD, 2MAD", sum(log_mcgf < mcg_3mads), sum(log_mcgf < mcg_2mads), "\n")
    hist(log_mcg, breaks=200, 
        xlab=expression(Log[10]~"MCG %"),
        ylab="No. cells", main=dirx)
    #abline(v = mcg_3mads, col='red', lwd=2)
    #abline(v = mcg_2mads, col='blue', lwd=2)
    #abline(v = mcg_1mads, col='purple', lwd=2)
    abline(v = mcg_med, col='darkgreen', lwd=2)
    abline(v = mcg_3madsp, col='red', lwd=2)
    abline(v = mcg_2madsp, col='blue', lwd=2)
    abline(v = mcg_1madsp, col='purple', lwd=2)
    dev.off()
    10^mcg_3madsp
}

plot_feat_counts_hist = function(per.feat, dfx, dirx, fname){
    featcts = per.feat$detected * dim(dfx)[2]/100.0
    log_featcts = log10(featcts)
    png(file=fname)
    xhist = hist(log_featcts, breaks=120,
            xlab=expression(Log[10]~"No. cells"), ylab="No. genes",
            main=dirx, prob=TRUE)
    lines(density(log_featcts), col="blue", lwd=2)
    xrange = (xhist$mids < 2) & (xhist$mids > 0)  & (xhist$counts > 0)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    abline(v = xrnmidv, col='red', lwd=2)
    log_fcts_filter = log_featcts > xrnmidv
    # cat("Feat : xrange, xrange_min, filter ", sum(xrange), xrnmidv, sum(log_featcts > xrnmidv), "\n")
    dev.off()
     
    10^(xrnmidv)
}

plot_feat_avg_hist = function(per.feat, dirx, fname){
    featavg = per.feat$mean 
    log_featavg = log10(featavg)
    png(file=fname)
    xhist = hist(log_featavg, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"), ylab="Prob. No. cells", 
        main=dirx, prob=TRUE)
    lines(density(log_featavg), col="blue", lwd=2)
    xrange = (xhist$mids < -1) & (xhist$mids > -4)  & (xhist$density > 0.05)
    xrngcts =xhist$counts[xrange]
    xrngmids = xhist$mids[xrange]
    xrnmidv = xrngmids[which.min(xrngcts)]
    abline(v = xrnmidv, col='red')
    #print(xrngcts)
    #print(xhist)
    log_favg_filter = log_featavg > xrnmidv
    # cat("Feat : xrange, xrange_min, filter ", sum(xrange), xrnmidv, sum(log_favg_filter), "\n")
    dev.off()
}

plot_counts_avg_hist = function(per.cell, dfx, dirx, fname){
    libmean = per.cell$sum / dim(dfx)[2]
    png(file=fname)
    libmeanf= libmean
    libmean = libmean[libmean <= 10]
    lmean_mad = mad(libmean)
    lmean_med = median(libmean)
    lmean_mean = mean(libmean)
    lmean_dev = sd(libmean)
    lmean_3mads = lmean_med - 3 * lmean_mad
    lmean_2mads = lmean_med - 2 * lmean_mad
    lmean_1mads = lmean_med - lmean_mad
    lmean_3madsp = lmean_med + 3 * lmean_mad
    lmean_2madsp = lmean_med + 2 * lmean_mad
    lmean_1madsp = lmean_med + lmean_mad
    lmean_3sds = lmean_mean - 3 * lmean_dev
    lmean_2sds = lmean_mean - 2 * lmean_dev
    lmean_1sds = lmean_mean - lmean_dev
    lmean_3sdsp = lmean_mean + 3 * lmean_dev
    lmean_2sdsp = lmean_mean + 2 * lmean_dev
    lmean_1sdsp = lmean_mean + lmean_dev
    #cat("median, mad, 3MAD, 2MAD", lmean_med, lmean_mad, sum(libmeanf < lmean_3mads), sum(libmeanf < lmean_2mads), "\n")
    hist(libmean, breaks=120, xlab=expression("Avg. No. Reads < 10"), ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = lmean_3mads, col='red', lwd=2)
    abline(v = lmean_2mads, col='blue', lwd=2)
    abline(v = lmean_1mads, col='purple', lwd=2)
    abline(v = lmean_med, col='darkgreen', lwd=2)
    abline(v = lmean_3madsp, col='red', lwd=2)
    abline(v = lmean_2madsp, col='blue', lwd=2)
    abline(v = lmean_1madsp, col='purple', lwd=2)
    
    abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    abline(v = lmean_med, col='darkgreen', lwd=2)
    dev.off()
}

plot_counts_logavg_hist = function(per.cell, dfx, dirx, fname){
    libmean = per.cell$sum / dim(dfx)[2]
    log_libmean = log10(libmean)
    png(file=fname)
    log_libmeanf= log_libmean
    #log_libmean = log_libmean[log_libmean < 1]
    lmean_mad = mad(log_libmean)
    lmean_med = median(log_libmean)
    lmean_mean = mean(log_libmean)
    lmean_dev = sd(log_libmean)
    lmean_3mads = lmean_med - 3 * lmean_mad
    lmean_2mads = lmean_med - 2 * lmean_mad
    lmean_1mads = lmean_med - lmean_mad
    lmean_3madsp = lmean_med + 3 * lmean_mad
    lmean_2madsp = lmean_med + 2 * lmean_mad
    lmean_1madsp = lmean_med + lmean_mad
    lmean_3sds = lmean_mean - 3 * lmean_dev
    lmean_2sds = lmean_mean - 2 * lmean_dev
    lmean_1sds = lmean_mean - lmean_dev
    lmean_3sdsp = lmean_mean + 3 * lmean_dev
    lmean_2sdsp = lmean_mean + 2 * lmean_dev
    lmean_1sdsp = lmean_mean + lmean_dev
    nlength = length(log_libmeanf)
    cat("Log Avg counts median, mad, 3MAD, 2MAD", 
      lmean_med, lmean_mad, 
      sum(log_libmeanf < lmean_3madsp), sum(log_libmeanf > lmean_1madsp),
      sum(log_libmeanf < lmean_3madsp & log_libmeanf > lmean_1madsp),
      sum(log_libmeanf < lmean_3madsp)/nlength, 
      sum(log_libmeanf > lmean_1madsp)/nlength,
      sum(log_libmeanf < lmean_3madsp & log_libmeanf > lmean_1madsp)/nlength,
       "\n")
    hist(log_libmean, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"), ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = lmean_3mads, col='red', lwd=2)
    abline(v = lmean_2mads, col='blue', lwd=2)
    abline(v = lmean_1mads, col='purple', lwd=2)
    abline(v = lmean_med, col='darkgreen', lwd=2)
    abline(v = lmean_3madsp, col='red', lwd=2)
    abline(v = lmean_2madsp, col='blue', lwd=2)
    abline(v = lmean_1madsp, col='purple', lwd=2)
    
    abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    abline(v = lmean_2sdsp, col='blue', lwd=2, lty=2)
    abline(v = lmean_1sdsp, col='purple', lwd=2, lty=2)
    dev.off()
    lmean_3madsp
}

plot_counts_logavg_hist2 = function(per.cell, lmean_3madsp, dfx, dirx, fname){
    libmean = per.cell$sum / dim(dfx)[2]
    log_libmean = log10(libmean)
    png(file=fname)
    log_libmeanf= log_libmean
    log_libmean = log_libmean[log_libmean < lmean_3madsp]
    lmean_mad = mad(log_libmean)
    lmean_med = median(log_libmean)
    lmean_mean = mean(log_libmean)
    lmean_dev = sd(log_libmean)
    lmean_3mads = lmean_med - 3 * lmean_mad
    lmean_2mads = lmean_med - 2 * lmean_mad
    lmean_1mads = lmean_med - lmean_mad
    lmean_3madsp = lmean_med + 3 * lmean_mad
    lmean_2madsp = lmean_med + 2 * lmean_mad
    lmean_1madsp = lmean_med + lmean_mad
    lmean_3sds = lmean_mean - 3 * lmean_dev
    lmean_2sds = lmean_mean - 2 * lmean_dev
    lmean_1sds = lmean_mean - lmean_dev
    lmean_3sdsp = lmean_mean + 3 * lmean_dev
    lmean_2sdsp = lmean_mean + 2 * lmean_dev
    lmean_1sdsp = lmean_mean + lmean_dev
    hist(log_libmean, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"),
	 ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = lmean_3mads, col='red', lwd=2)
    abline(v = lmean_2mads, col='blue', lwd=2)
    abline(v = lmean_1mads, col='purple', lwd=2)
    abline(v = lmean_med, col='darkgreen', lwd=2)
    abline(v = lmean_3madsp, col='red', lwd=2)
    abline(v = lmean_2madsp, col='blue', lwd=2)
    abline(v = lmean_1madsp, col='purple', lwd=2)
    
    abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    abline(v = lmean_2sdsp, col='blue', lwd=2, lty=2)
    abline(v = lmean_1sdsp, col='purple', lwd=2, lty=2)
    dev.off()
}

plot_ngenes = function(per.cell, dirx, fname){
    ngenes = per.cell$detected
    log_ngenes = log10(ngenes)
    png(file=fname)
    log_ngenesf = log_ngenes
    #log_ngenes = log_ngenes[log_ngenes < 5]
    ngenes_nmad = isOutlier(log_ngenes, nmads=3, typ="higher")
    ngenes_mad = mad(log_ngenes)
    ngenes_med = median(log_ngenes)
    ngenes_3mads = ngenes_med - 3 * ngenes_mad
    ngenes_2mads = ngenes_med - 2 * ngenes_mad
    ngenes_1mads = ngenes_med -  ngenes_mad
    ngenes_3madsp = ngenes_med + 3 * ngenes_mad
    ngenes_2madsp = ngenes_med + 2 * ngenes_mad
    ngenes_1madsp = ngenes_med +  ngenes_mad
    nlength = length(log_ngenesf)
    cat("Log Ngenes : median, mad, <3MADP, >1MAD Both",
      ngenes_med, ngenes_mad, 
      sum(log_ngenesf < ngenes_3madsp), 
      sum(log_ngenesf > ngenes_1mads),
      sum(log_ngenesf < ngenes_3madsp & log_ngenesf > ngenes_1mads),
      sum(log_ngenesf < ngenes_3madsp)/nlength,
      sum(log_ngenesf > ngenes_1mads)/nlength,
      sum(log_ngenesf < ngenes_3madsp & log_ngenesf > ngenes_1mads)/nlength
        "\n")
    hist(log_ngenesf, breaks=120, xlab=expression(Log[10]~"No. genes"), ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = ngenes_3mads, col='red', lwd=2)
    abline(v = ngenes_2mads, col='blue', lwd=2)
    abline(v = ngenes_1mads, col='purple', lwd=2)
    abline(v = ngenes_med, col='darkgreen', lwd=2)
    abline(v = ngenes_3madsp, col='red', lwd=2)
    abline(v = ngenes_2madsp, col='blue', lwd=2)
    abline(v = ngenes_1madsp, col='purple', lwd=2)
    dev.off()
}

plot_hist = function(per.cell, per.feat, dfx, dirx, out_dir, out_prefix){
    cat("Dims: ",  dirx, dim(dfx), "\n")

    fname = paste(out_dir, dirx, 
        paste(out_prefix,
            "-scater-qc-mcgpct-hist.png", sep=""),
        sep="/"  )
    mcg_threshold = plot_mcg_hist(per.cell, dirx, fname)

    fname = paste(out_dir, dirx, 
        paste(out_prefix, 
            "-scater-qc-feat-hist.png", sep=""),
        sep="/"  )
    feat_cells_threshold = plot_feat_counts_hist(per.feat, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-feats-avg-hist.png", sep=""),
        sep="/"  )
    plot_feat_avg_hist(per.feat, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-counts-avg-hist.png", sep=""),
        sep="/"  )
    plot_counts_avg_hist(per.cell, dfx, dirx, fname)

    fname = paste(out_dir, dirx, 
        paste(out_prefix, "-scater-qc-counts-logavg-hist.png", sep=""),
        sep="/"  )
    lmean_3madsp = plot_counts_logavg_hist(per.cell, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-counts-logavg-hist2.png", sep=""),
        sep="/"  )
    plot_counts_logavg_hist2(per.cell, lmean_3madsp, dfx, dirx, fname)

    fname = paste(out_dir, dirx,
        paste(out_prefix, "-scater-qc-ngenes-hist.png", sep=""),
        sep="/"  )
    plot_ngenes(per.cell, dirx, fname)

    c(mcg=mcg_threshold, feature=feat_cells_threshold)
}

plot_qcstats = function(out_dir, out_prefix, in_dirs){
    for(dirx in in_dirs){
        dfx = read10xCounts(dirx)
        per.cell = perCellQCMetrics(dfx, 
            subset=list(MCG=grep("AT[MC]G",
                            rownames(dfx))))
        per.feat = perFeatureQCMetrics(dfx)
        plot_flip_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        # all_thrs = plot_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        # seurat_qcplot(dirx, all_thrs, out_dir, out_prefix, ".png")
    }
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) >= 3){
    plot_qcstats(args[1], args[2], args[3:length(args)])
}  else {
    print(args)
    print("Usage: Rscript scater_qcplot.R outdir out_prefix indir1 indir2 ...")
}
