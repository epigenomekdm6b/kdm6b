#!/usr/bin/Rscript

args <- commandArgs(TRUE)
if(length(args) < 4) {
        cat("\nUsage: r_seurat_extract_then_subclustering_by_cca.R [PROJECT_HOME_ID] [DIR_SO_PREVIOUS] [SO_PREFIX] [TARGET_CLUSTERS]
        ./r_seurat_extract_then_subclustering.R mll4_e15_arcuate_nucleus out_mll4_e15_arcuate_nucleus_rem_hb clu0k7 0,4\n\n")
        q()
}

dims.use <- 5

prj.home.id <- args[1]
dir.so.prev <- args[2]
so.prefix <- args[3]
str.target.clusters <- args[4]

prj.id <- paste(prj.home.id,so.prefix,sep="_")

dir.out <- paste('out',prj.home.id,so.prefix,sep='_')
system(paste('rm -rf',dir.out))
system(paste('mkdir',dir.out))

dir.before <- getwd() #print(dir.before)

setwd(dir.out)

# Load the dataset

path.so <- paste('..',dir.so.prev,'so.rdb',sep='/')

# temporary!!!
load(path.so)
#print(path.so)
#q()
library(Seurat)
library(dplyr)

#----------------------------
# Extract barcodes from clusters

target.clusters <- unlist(strsplit(str.target.clusters,','))
#target.clusters <- c( "0")

extract.barcodes <- c()
for(i in 1:length(target.clusters)) {
        target.cluster <- target.clusters[i]
        target.barcodes <- names(so@ident[so@ident == target.cluster])
        extract.barcodes <- c(extract.barcodes, target.barcodes)
}

length(extract.barcodes)

##################################################

subset.matrix <- so@raw.data[,extract.barcodes]

# create two independant 'so' by label 's1' and 's2'

ix.s1 <- grep("s1_",colnames(subset.matrix))
ix.s2 <- grep("s2_",colnames(subset.matrix))

subset.matrix.s1 <- subset.matrix[,ix.s1]
subset.matrix.s2 <- subset.matrix[,ix.s2]

so.s1 <- CreateSeuratObject(raw.data = subset.matrix.s1, min.cells = 3, min.genes = 200, project = 's1')
so.s2 <- CreateSeuratObject(raw.data = subset.matrix.s2, min.cells = 3, min.genes = 200, project = 's2')

#-----------------------------------------------------------------------------
mito.genes <- grep(pattern = "^mt-", x = rownames(x = so.s1@data), value = TRUE)
percent.mito <- Matrix::colSums(so.s1@raw.data[mito.genes, ])/Matrix::colSums(so.s1@raw.data)
so.s1 <- AddMetaData(object = so.s1, metadata = percent.mito, col.name = "percent.mito")
#-----------------------------------------------------------------------------
mito.genes <- grep(pattern = "^mt-", x = rownames(x = so.s2@data), value = TRUE)
percent.mito <- Matrix::colSums(so.s2@raw.data[mito.genes, ])/Matrix::colSums(so.s2@raw.data)
so.s2 <- AddMetaData(object = so.s2, metadata = percent.mito, col.name = "percent.mito")
#-----------------------------------------------------------------------------

###############################
#low.thresholds <- 200
#high.thresholds <- 6000
###############################
#so.s1 <- FilterCells(object = so.s1, subset.names = c("nGene", "percent.mito"),
    #low.thresholds = c(low.thresholds, -Inf), high.thresholds = c(high.thresholds, 0.05))

so.s1 <- NormalizeData(object = so.s1, normalization.method = "LogNormalize", scale.factor = 10000) 
png('tmp.png')
so.s1 <- FindVariableGenes(object = so.s1, mean.function = ExpMean, dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()

length(x = so.s1@var.genes)
#[1] 1891

so.s1 <- ScaleData(object = so.s1, vars.to.regress = c("nUMI", "percent.mito"))
#-----------------------------------------------------------------------------
#so.s2 <- FilterCells(object = so.s2, subset.names = c("nGene", "percent.mito"),
    #low.thresholds = c(low.thresholds, -Inf), high.thresholds = c(high.thresholds, 0.05))

so.s2 <- NormalizeData(object = so.s2, normalization.method = "LogNormalize", scale.factor = 10000) 
png('tmp2.png')
so.s2 <- FindVariableGenes(object = so.s2, mean.function = ExpMean, dispersion.function = LogVMR,
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = so.s2@var.genes)
#[1] 1891

so.s2 <- ScaleData(object = so.s2, vars.to.regress = c("nUMI", "percent.mito"))
#-----------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# common in 'cca'

cl.s1 <- 'CTRL'
cl.s2 <- 'KO'

# Set up control object
so.s1@meta.data$stim <- cl.s1
so.s2@meta.data$stim <- cl.s2

# Gene selection for input to CCA
png('tmp2.png')
so.s1 <- FindVariableGenes(so.s1, do.plot = F)
dev.off()
png('tmp2.png')
so.s2 <- FindVariableGenes(so.s2, do.plot = F)
dev.off()

g.1 <- so.s1@var.genes
g.2 <- so.s2@var.genes
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(so.s1@scale.data))
genes.use <- intersect(genes.use, rownames(so.s2@scale.data))

#-----------------------------
# Perform a canonical correlation analysis (CCA)
#-----------------------------

combined <- RunCCA(so.s1, so.s2, genes.use = genes.use, num.cc = 30)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "stim",
    pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "stim",
    do.return = TRUE)

file.png <- paste(so.prefix, "_plot_1_pre.png", sep="")
png(file.png) 
plot_grid(p1, p2)
dev.off()

PrintDim(object = combined, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

#-----------------------------
# Align the CCA subspaces
#-----------------------------

combined <- AlignSubspace(combined, reduction.type = "cca", grouping.var = "stim",
    dims.align = 1:dims.use)

p1 <- VlnPlot(object = combined, features.plot = "ACC1", group.by = "stim",
    do.return = TRUE)
p2 <- VlnPlot(object = combined, features.plot = "ACC2", group.by = "stim",
    do.return = TRUE)

file.png <- paste(so.prefix, "_plot_4_pre.png", sep="")
png(file.png) 
plot_grid(p1, p2)
dev.off()

#-----------------------------
# Perform an integrated analysis
#-----------------------------

# t-SNE and Clustering
combined <- RunTSNE(combined, reduction.use = "cca.aligned", dims.use = 1:dims.use,
    do.fast = T)
combined <- FindClusters(combined, reduction.type = "cca.aligned",
    resolution = 0.6, dims.use = 1:dims.use)

# Visualization
p1 <- TSNEPlot(combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(combined, do.label = T, do.return = T, pt.size = 0.5)

file.png <- paste(so.prefix, "_plot_5_tsne.png", sep="")
png(file.png) 
plot_grid(p1, p2)
dev.off()

file.pdf <- paste(so.prefix, "_plot_5_tsne.pdf", sep="")
pdf(file=file.pdf, width=14, height=7)
plot_grid(p1, p2)
dev.off()

#************************
so <- combined
#************************
#---------------------------------------------
so.table <- table(so@ident)

file.cluster.number <- paste(so.prefix, '_cluster_number.txt',sep="")
write.table(file=file.cluster.number,so.table,quote=F,row.names=F,col.names=F,sep="\t")

#---------------------------------------------
th.min.pct <- 0
so.markers <- FindAllMarkers(object = so, only.pos = FALSE, min.pct = th.min.pct)

file.markers_cluster <- paste(so.prefix, '_markers_cluster.tab', sep="")
write.table(file=file.markers_cluster,so.markers,sep="\t",quote=F)

#---------------------------------------------
#abby
which.barcodes <- function(seurat.object, target.cluster) {
        target.barcodes <- names(seurat.object@ident[seurat.object@ident == target.cluster])
        return(target.barcodes)
}

genes.in.order <- row.names(so@raw.data)

clusters.all <- levels(so@ident)

n.clusters <- length(clusters.all)
n.genes <- length(genes.in.order)

mtx <- matrix(NA,n.genes,n.clusters)

for(i in 1:n.clusters) {

        cluster.target <- clusters.all[i]
        barcodes.target <- which.barcodes(so, cluster.target)
        mean.cluster <- apply(so@data[,barcodes.target], 1, mean)
        mtx[,i] <- mean.cluster
}

df.mean <- as.data.frame(mtx)
row.names(df.mean) <- row.names(so@raw.data)
names(df.mean) <- clusters.all
file.mean.expression <- paste(so.prefix,"_mean_expression_in_cluster_data.tab",sep="")
write.table(file=file.mean.expression,df.mean,sep="\t",quote=F)

save(so,file="so.rdb")







































