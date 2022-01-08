#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if(length(args) < 1) {
        cat("\nUsage: r_seurat_split_cell_from_so_by_condition_from_pca.R [SO_PREFIX] 
        ./r_seurat_split_cell_from_so_by_condition_from_pca.R ctr\n\n")
        q()
}

so.prefix <- args[1]

library(Seurat)
load('so.rdb')

so <- RunPCA(object = so, pc.genes = so@var.genes, do.print = TRUE, pcs.print = 1:5,
    genes.print = 10)

# Examine and visualize PCA results a few different ways
PrintPCA(object = so, pcs.print = 1:10, genes.print = 5, use.full = FALSE)

#----------------------------------------------------
file.png <- paste(so.prefix, "_plot_pca.png", sep="")
png(file.png)
PCAPlot(object = so, dim.1 = 1, dim.2 = 2)
dev.off()

file.pdf <- paste(so.prefix, "_plot_pca.pdf", sep="")
pdf(file.pdf)
PCAPlot(object = so, dim.1 = 1, dim.2 = 2)
dev.off()

#----------------------------------------------------

so <- FindClusters(object = so, reduction.type = "pca", dims.use = 1:10,
    resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = so)

so.table <- table(so@ident)

file.cluster.number <- paste(so.prefix, '_cluster_number.txt',sep="")
write.table(file=file.cluster.number,so.table,quote=F,row.names=F,col.names=F,sep="\t")

so <- RunTSNE(object = so, dims.use = 1:10, do.fast = TRUE)

#-----------------------
file.png <- paste(so.prefix, "_plot_tsne.png", sep="")
png(file.png)
#---------------------
TSNEPlot(object = so,do.label=T)
#---------------------
dev.off()

file.pdf <- paste(so.prefix, "_plot_tsne.pdf", sep="")
pdf(file.pdf)
#---------------------
TSNEPlot(object = so,do.label=T)
#---------------------
dev.off()

th.min.pct <- 0
so.markers <- FindAllMarkers(object = so, only.pos = FALSE, min.pct = th.min.pct)
#so.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

file.markers_cluster <- paste(so.prefix, '_markers_cluster.tab', sep="")
write.table(file=file.markers_cluster,so.markers,sep="\t",quote=F)


#---------------------------------------------------------------
#abby
which.barcodes <- function(seurat.object, target.cluster) {
        target.barcodes <- names(seurat.object@ident[seurat.object@ident == target.cluster])
        return(target.barcodes)
}

genes.in.order <- row.names(so@raw.data)
########################################################

clusters.all <- levels(so@ident)

n.clusters <- length(clusters.all)
n.genes <- length(genes.in.order)

mtx <- matrix(NA,n.genes,n.clusters)

for(i in 1:n.clusters) {

        cluster.target <- clusters.all[i]
        barcodes.target <- which.barcodes(so, cluster.target)
        ########################################################
        mean.cluster <- apply(so@data[,barcodes.target], 1, mean)
        ########################################################
        mtx[,i] <- mean.cluster
}

df.mean <- as.data.frame(mtx)
########################################################
row.names(df.mean) <- row.names(so@raw.data)
names(df.mean) <- clusters.all
########################################################

########################################################
#file.mean.expression <- paste(so.prefix,"_mean_expression_in_cluster_scale_data.tab",sep="")
file.mean.expression <- paste(so.prefix,"_mean_expression_in_cluster_data.tab",sep="")
########################################################
write.table(file=file.mean.expression,df.mean,sep="\t",quote=F)
#---------------------------------------------------------------

save(so,so.markers, file="so.rdb")

