library(Matrix)
gen_subset= function(matrix_file, cols_file, rows_file, filter_col_file, output_file){
    mx = readMM(matrix_file)
    cxdf = read.table(cols_file, stringsAsFactors=F, header=F)
    rxdf = read.table(rows_file, stringsAsFactors=F, header=F)
    mdf = as.data.frame(as.matrix(mx))
    rownames(mdf) = rxdf$V1
    colnames(mdf) = cxdf$V1
    fndf = read.table(filter_col_file, stringsAsFactors=F, header=F)

    fxdf = mdf[,  colnames(mdf) %in% fndf$V1]
    write.csv(fxdf, output_file)
}

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 5){
    gen_subset(args[1], args[2], args[3], args[4], args[5])
} 
if(length(args) != 5){
    cat("Required args: matrix_file, cols_file, rows_file, filter_col_file, output_file")
}