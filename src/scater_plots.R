library(scater)
library("DropletUtils")

qcstat_dir = function(dirx){
    dfx = read10xCounts(dirx)
    per.cell = perCellQCMetrics(dfx, 
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    per.feat = perFeatureQCMetrics(dfx)
    tfx = per.feat$detected * dim(dfx)[2]/100.0
    libsize_mads = isOutlier(per.cell$sum, nmads=3,
                            type="lower", log=FALSE)
    genects_mads = isOutlier(per.cell$detected, nmads=3,
                            type="lower", log=FALSE)
    feat_cell_cts = per.feat$detected * dim(dfx)[2]/100.0
    feat_cell_mean = per.feat$mean >= 1.0

    libsize_95p = per.cell$sum < (quantile(per.cell$sum, .95) * 0.05)
    genects_99p = per.cell$detected < (quantile(per.cell$detected, 0.99) * 0.05)

    mcg_nmads = isOutlier(per.cell$subsets_MCG_percent, nmads=3, typ="higher")
    mcg_5p = per.cell$subsets_MCG_percent > 5
    mcg_4p = per.cell$subsets_MCG_percent > 4
    mcg_3p = per.cell$subsets_MCG_percent > 3
    c(#dirname = dirx,
      total_cells = dim(dfx)[2],
      libsize_mads=sum(libsize_mads),
      genects_mads=sum(genects_mads),
      libsize_95p=sum(libsize_95p),
      genects_99p=sum(genects_99p),
      mcg_nmads=sum(mcg_nmads),
      mcg_5p=sum(mcg_5p),
      mcg_4p=sum(mcg_4p),
      mcg_3p=sum(mcg_3p),
      total_features = dim(dfx)[1], 
      feat_cell_mean=sum(feat_cell_mean)
       )

}

get_qcstats = function(outfile, indirs){
    print(outfile)
    print(indirs)
    qcl = lapply(indirs, qcstat_dir)
    qdf = t(as.data.frame(qcl))
    rownames(qdf) = indirs
    write.table(qdf, outfile, sep="\t", quote=FALSE)
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

plot_feat_avg_hist = function(per.feat, dirx, fname, midv=TRUE){
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
    if (midv == TRUE){ abline(v = xrnmidv, col='red') }
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
    
    #abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    #abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    #abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    #abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    #abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    #abline(v = lmean_med, col='darkgreen', lwd=2)
    dev.off()
}

plot_logcounts_hist = function(per.cell, dfx, dirx, fname, print_values=TRUE){
    libsize = per.cell$sum 
    log_libsize = log10(libsize)
    png(file=fname)
    log_libsizef= log_libsize
    #log_libsize = log_libsize[log_libsize < 1]
    lsize_mad = mad(log_libsize)
    lsize_med = median(log_libsize)
    lsize_mean = mean(log_libsize)
    lsize_dev = sd(log_libsize)
    lsize_3mads = lsize_med - 3 * lsize_mad
    lsize_2mads = lsize_med - 2 * lsize_mad
    lsize_1mads = lsize_med - lsize_mad
    lsize_3madsp = lsize_med + 3 * lsize_mad
    lsize_2madsp = lsize_med + 2 * lsize_mad
    lsize_1madsp = lsize_med + lsize_mad
    lsize_3sds = lsize_mean - 3 * lsize_dev
    lsize_2sds = lsize_mean - 2 * lsize_dev
    lsize_1sds = lsize_mean - lsize_dev
    lsize_3sdsp = lsize_mean + 3 * lsize_dev
    lsize_2sdsp = lsize_mean + 2 * lsize_dev
    lsize_1sdsp = lsize_mean + lsize_dev
    nlength = length(log_libsizef)
    if (print_values == TRUE){
     cat("Log Avg counts ", 
      "median : ", 10^lsize_med, 10^lsize_mad,
      "LB avg cts: ", 10^lsize_1mads, 
      "UB avg cts: ", 10^lsize_3madsp,
      "LB ncells:", sum(log_libsizef > lsize_1mads),
      sum(log_libsizef > lsize_1mads)*100/nlength,
      "UB ncells:", sum(log_libsizef < lsize_3madsp), 
      sum(log_libsizef < lsize_3madsp)*100/nlength, 
      "UB+LB ncells:", sum(log_libsizef < lsize_3madsp & log_libsizef > lsize_1mads),
      sum(log_libsizef < lsize_3madsp & log_libsizef > lsize_1mads)*100/nlength,
       "\n")
    }
    hist(log_libsize, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"), ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = lsize_3mads, col='red', lwd=2)
    abline(v = lsize_2mads, col='blue', lwd=2)
    abline(v = lsize_1mads, col='purple', lwd=2)
    abline(v = lsize_med, col='darkgreen', lwd=2)
    abline(v = lsize_3madsp, col='red', lwd=2)
    abline(v = lsize_2madsp, col='blue', lwd=2)
    abline(v = lsize_1madsp, col='purple', lwd=2)
    
    #abline(v = lsize_3sds, col='red', lwd=2, lty=2)
    #abline(v = lsize_2sds, col='blue', lwd=2, lty=2)
    #abline(v = lsize_1sds, col='purple', lwd=2, lty=2)
    #abline(v = lsize_mean, col='darkgreen', lwd=2, lty=2)
    #abline(v = lsize_3sdsp, col='red', lwd=2, lty=2)
    #abline(v = lsize_2sdsp, col='blue', lwd=2, lty=2)
    #abline(v = lsize_1sdsp, col='purple', lwd=2, lty=2)
    dev.off()
    lsize_3madsp
}



plot_counts_logavg_hist = function(per.cell, dfx, dirx, fname, print_values=TRUE){
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
    if (print_values == TRUE){
     cat("Log Avg counts ", 
      "median : ", 10^lmean_med, 10^lmean_mad,
      "LB avg cts: ", 10^lmean_1mads, 
      "UB avg cts: ", 10^lmean_3madsp,
      "LB ncells:", sum(log_libmeanf > lmean_1mads),
      sum(log_libmeanf > lmean_1mads)*100/nlength,
      "UB ncells:", sum(log_libmeanf < lmean_3madsp), 
      sum(log_libmeanf < lmean_3madsp)*100/nlength, 
      "UB+LB ncells:", sum(log_libmeanf < lmean_3madsp & log_libmeanf > lmean_1mads),
      sum(log_libmeanf < lmean_3madsp & log_libmeanf > lmean_1mads)*100/nlength,
       "\n")
    }
    hist(log_libmean, breaks=120, xlab=expression(Log[10]~"Avg. No. Reads"), ylab="No. cells", main=dirx, prob=FALSE)
    abline(v = lmean_3mads, col='red', lwd=2)
    abline(v = lmean_2mads, col='blue', lwd=2)
    abline(v = lmean_1mads, col='purple', lwd=2)
    abline(v = lmean_med, col='darkgreen', lwd=2)
    abline(v = lmean_3madsp, col='red', lwd=2)
    abline(v = lmean_2madsp, col='blue', lwd=2)
    abline(v = lmean_1madsp, col='purple', lwd=2)
    
    #abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    #abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    #abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    #abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    #abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    #abline(v = lmean_2sdsp, col='blue', lwd=2, lty=2)
    #abline(v = lmean_1sdsp, col='purple', lwd=2, lty=2)
    dev.off()
    lmean_3madsp
}

plot_counts_logavg_hist2 = function(per.cell, lmean_3madsp, dfx, dirx, fname) {
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
    
    #abline(v = lmean_3sds, col='red', lwd=2, lty=2)
    #abline(v = lmean_2sds, col='blue', lwd=2, lty=2)
    #abline(v = lmean_1sds, col='purple', lwd=2, lty=2)
    #abline(v = lmean_mean, col='darkgreen', lwd=2, lty=2)
    #abline(v = lmean_3sdsp, col='red', lwd=2, lty=2)
    #abline(v = lmean_2sdsp, col='blue', lwd=2, lty=2)
    #abline(v = lmean_1sdsp, col='purple', lwd=2, lty=2)
    dev.off()
}

plot_ngenes = function(per.cell, dirx, fname, print_values=TRUE){
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
    if (print_values == TRUE){
     cat("Log Ngenes : ",
      "Median:", 10^ngenes_med, 10^ngenes_mad, 
      "lB ngenes : ", 10^ngenes_1mads, 
      "uB ngenes : ", 10^ngenes_3madsp,
      "LB no. cells", sum(log_ngenesf > ngenes_1mads),
      sum(log_ngenesf > ngenes_1mads)*100/nlength,
      "UB no. cells: ", sum(log_ngenesf < ngenes_3madsp), 
      sum(log_ngenesf < ngenes_3madsp)*100/nlength,
      "LB+UB no. cells: ", sum(log_ngenesf < ngenes_3madsp & log_ngenesf > ngenes_1mads),
      sum(log_ngenesf < ngenes_3madsp & log_ngenesf > ngenes_1mads)/nlength,
        "\n")
    }
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
    outfile = paste(out_dir, "scater-qc.tsv", sep="/")
    get_qcstats(outfile, indirs)
    for(dirx in in_dirs){
        dfx = read10xCounts(dirx)
        per.cell = perCellQCMetrics(dfx, 
            subset=list(MCG=grep("AT[MC]G",
                            rownames(dfx))))
        per.feat = perFeatureQCMetrics(dfx)
        #plot_flip_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        all_thrs = plot_hist(per.cell, per.feat, dfx, dirx, out_dir, out_prefix)
        #seurat_ba_qcplot(dirx, all_thrs, out_dir, out_prefix, ".png")
    }
}


plot_cells_hist = function(dfx, dirx, name_prefix,
                           out_dir, out_prefix){
    per.cell = perCellQCMetrics(dfx,
        subset=list(MCG=grep("AT[MC]G",
                        rownames(dfx))))
    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
              "-scater-qc-logcounts-hist.png", sep=""),
        sep="/"  )
    # print(fname)
    plot_logcounts_hist(per.cell, dfx, dirx, fname, FALSE)

    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
              "-scater-qc-counts-logavg-hist.png", sep=""),
        sep="/"  )
    # print(fname)
    plot_counts_logavg_hist(per.cell, dfx, dirx, fname, FALSE)

    fname = paste(out_dir, dirx,
        paste(out_prefix, name_prefix,
          "-scater-qc-ngenes-hist.png", sep=""),
        sep="/"  )
    plot_ngenes(per.cell, dirx, fname, FALSE)
}
