library(Seurat)
load('so.rdb')
so.prefix <- 'comp_rep2_control_vs_ko_slc18a3'

th.min.pct <- 0
so.markers <- FindAllMarkers(object = so, only.pos = FALSE, min.pct = th.min.pct)

file.markers_cluster <- paste(so.prefix, '_markers_cluster.tab', sep="")
write.table(file=file.markers_cluster,so.markers,sep="\t",quote=F)

#--------------------------------------------------------------------------------------------------

which.barcodes <- function(seurat.object, target.cluster) {
	        target.barcodes <- names(seurat.object@ident[seurat.object@ident == target.cluster])
        return(target.barcodes)
}

#genes.in.order <- row.names(so@raw.data)
genes.in.order <- row.names(so@data)

clusters.all <- levels(so@ident)

n.clusters <- length(clusters.all)
n.genes <- length(genes.in.order)

mtx <- matrix(NA,n.genes,n.clusters)

for(i in 1:n.clusters) {

	        cluster.target <- clusters.all[i]
        	barcodes.target <- which.barcodes(so, cluster.target)
	        #mean.cluster <- apply(so@raw.data[,barcodes.target], 1, mean)
	        mean.cluster <- apply(so@data[,barcodes.target], 1, mean)
	        mtx[,i] <- mean.cluster
}

df.mean <- as.data.frame(mtx)
#row.names(df.mean) <- row.names(so@raw.data)
row.names(df.mean) <- row.names(so@data)
names(df.mean) <- clusters.all

#file.mean.expression <- paste(so.prefix,"_mean_expression_in_cluster_raw_data.tab",sep="")
file.mean.expression <- paste(so.prefix,"_mean_expression_in_cluster_data.tab",sep="")
write.table(file=file.mean.expression,df.mean,sep="\t",quote=F)

#--------------------------------------------------------------------------------------------

so.table <- table(so@ident)

file.cluster.number <- paste(so.prefix, '_cluster_number.txt',sep="")
write.table(file=file.cluster.number,so.table,quote=F,row.names=F,col.names=F,sep="\t")


