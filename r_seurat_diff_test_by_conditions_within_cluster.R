#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if(length(args) < 1) {
        cat("\nUsage: r_seurat_diff_test_by_conditions_within_cluster.R [PREFIX]
        ./r_seurat_diff_test_by_conditions_within_cluster.R comp_rep1_control_vs_ko\n\n")
        q()
}

#----------------------------------------------------------
SetIfNull <- function (x, default) 
#----------------------------------------------------------
{
    if (is.null(x = x)) { return(default) }
    else { return(x) }
}

#----------------------------------------------------------
file.clusters <- function(so.prefix, data.whole, clu.target) {
#----------------------------------------------------------
#for test
clu.target <- '9'

	#cells.in.clu <- WhichCells(so, ident=clu.target) # for Seurat v2.x
	cells.in.clu <- WhichCells(so, idents=clu.target) # for Seurat v4.x

	cell.names.with.label <- function(cell.names, prefix.label) {
        	return(cell.names[grep(prefix.label,cell.names)])
	}

	#cells.ctrl <- cell.names.with.label(cells.in.clu, "s1")   # for rep2
	#cells.case <- cell.names.with.label(cells.in.clu, "s2")   # for rep2
	cells.ctrl <- cell.names.with.label(cells.in.clu, "_1$")   # for rep1
	cells.case <- cell.names.with.label(cells.in.clu, "_2$")   # for rep1

	#==================================================
	res.diff.test <- WilcoxDETest(object = so,
            	cells.1 = cells.ctrl, cells.2 = cells.case, genes.use = genes.use,
            	print.bar = T)

	p.val <- res.diff.test$p_val
	#==================================================
	df.to.return <- data.frame()
	pseudocount.use=1

	#-----------------------------------------
	# mean expression & logFC
    	mean.data.1 <- apply(X = data.whole[genes.use, cells.ctrl, drop = F],
        	MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
            	pseudocount.use))
    	mean.data.2 <- apply(X = data.whole[genes.use, cells.case, drop = F],
        	MARGIN = 1, FUN = function(x) log(x = mean(x = expm1(x = x)) +
            	pseudocount.use))

        mean.ctrl <- expm1(mean.data.1)
        mean.case <- expm1(mean.data.2)
        log2FC <- log2(mean.case/mean.ctrl)
        log2FC[which(log2FC==-Inf)] = NA
        log2FC[which(log2FC==Inf)] = NA

    	#df.to.return$p_val_adj = p.adjust(p = df.to.return$p_val, method = "fdr",
    	p.adj = p.adjust(p = p.val, method = "fdr",
        	n = nrow(x = GetAssayData(object = so,
            	slot = "data")))
	#-----------------------------------------

	#-----------------------------------------
	# percentage of cells
    	thresh.min <- 0
    	pct.data.1 <- round(x = apply(X = data.whole[genes.use, cells.ctrl,
        	drop = F], MARGIN = 1, FUN = function(x) {
        	return(sum(x > thresh.min)/length(x = x))
    	}), digits = 3)
    	pct.data.2 <- round(x = apply(X = data.whole[genes.use, cells.case,
        	drop = F], MARGIN = 1, FUN = function(x) {
        	return(sum(x > thresh.min)/length(x = x))
    	}), digits = 3)
	#-----------------------------------------

	df.to.return <- data.frame(
		mean.ctrl=mean.ctrl,
		mean.case=mean.case,
		log2FC=log2FC,
		p.val=p.val,
		p.adj=p.adj,
		pct.ctrl=pct.data.1,
		pct.case=pct.data.2
	)

	out.file <- paste0(so.prefix,'_within_clu_',clu.target,'.tab')

	write.table(file=out.file,df.to.return,quote=F,sep="\t")

}
#----------------------------------------------------------

load('so.rdb')
library(Seurat)

#so.prefix <- args[1]
#so.prefix <- 'comp_rep2_control_vs_ko'
so.prefix <- 'comb_rep1_slc18a3'

dir.out <- paste('diff_test',so.prefix,sep='_')
system(paste('rm -rf',dir.out))
system(paste('mkdir',dir.out))
setwd(dir.out)

data.whole <- GetAssayData(object=so)
genes.use <- rownames(data.whole)

#----------------
#so.prefix <- 'comp_rep2_control_vs_ko'
#clu.target <- '1'
#----------------

for( i in 1:length(levels(so@ident))) {
	clu.target <- levels(so@ident)[i]
	file.clusters(so.prefix, data.whole, clu.target)
}
