#!/usr/bin/Rscript

# common until here...
##########################################################################################
# [ Cell number distribution per cluster ]

args <- commandArgs(TRUE)
if(length(args) < 1) {
        cat("\nUsage: r_seurat_overall_cell_distribution.R [SO_PREFIX]
        ./r_seurat_overall_cell_distribution.R clu_10\n\n")
        q()
}

so.prefix <- args[1]

options(digits=2)

out.file <- paste0(so.prefix, "_cell_distrib_per_cond.tab")

library(Seurat)

load('so.rdb')

clu.all <- levels(so@ident)
n.clus <- length(clu.all)

n.s1s <- c()
n.s2s <- c()
totals <- c()

r.s1s <- c()
r.s2s <- c()
for(i in 1:n.clus) { 
	
	cl.target <- clu.all[i]

	n.s1 <- length(grep("s1_",names(so@ident[so@ident==cl.target])))
	n.s2 <- length(grep("s2_",names(so@ident[so@ident==cl.target])))
	total <- n.s1 + n.s2

	r.s1 <- round(n.s1/total, digits=3)
	r.s2 <- round(n.s2/total, digits=3)

	n.s1s[i] <- n.s1	
	n.s2s[i] <- n.s2	
	totals[i] <- total
	r.s1s[i] <- r.s1	
	r.s2s[i] <- r.s2	

}

df.stat <- data.frame(Cluster=clu.all, N.s1=n.s1s, N.s2=n.s2s, Total=totals, R.s1=r.s1s, R.s2=r.s2s)
write.table(df.stat, file=out.file, quote=F,row.names=F,sep="\t")
