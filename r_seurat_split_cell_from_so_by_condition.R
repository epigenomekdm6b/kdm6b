#!/usr/bin/Rscript

#args <- commandArgs(TRUE)
#if(length(args) < 4) {
	        #cat("\nUsage: r_seurat_extract_then_subclustering_by_cca.R [PROJECT_HOME_ID] [DIR_SO_PREVIOUS] [SO_PREFIX] [TARGET_CLUSTERS]
		            #./r_seurat_extract_then_subclustering.R mll4_e15_arcuate_nucleus out_mll4_e15_arcuate_nucleus_rem_hb clu0k7 0,4\n\n")
			            #q()
#}

#prj.home.id <- args[1]
#dir.so.prev <- args[2]
#so.prefix <- args[3]
#str.target.clusters <- args[4]

prj.home.id <- 'comp_rep2_control_vs_ko_slc18a3'
dir.so.prev <- 'comp_rep2_control_vs_ko_slc18a3' 
so.prefix <- 'splitted'
#str.target.clusters <- '10'

prj.id <- paste(prj.home.id,so.prefix,sep="_")

#print(prj.id); q()

dir.out <- paste('out',prj.home.id,so.prefix,sep='_')
system(paste('rm -rf',dir.out))
system(paste('mkdir',dir.out))

dir.before <- getwd() #print(dir.before)

setwd(dir.out)
#===========================================================================================

# Load the dataset

path.so <- paste('..',dir.so.prev,'so.rdb',sep='/')

# temporary!!!
load(path.so)
#print(path.so)
#q()
library(Seurat)

library(dplyr)

#----------------------------
ix.s1 <- grep("s1_",so@cell.names)
ix.s2 <- grep("s2_",so@cell.names)

subset.matrix.s1 <- so@raw.data[,ix.s1]
subset.matrix.s2 <- so@raw.data[,ix.s2]

so.s1 <- CreateSeuratObject(raw.data = subset.matrix.s1, project = 's1')
so.s2 <- CreateSeuratObject(raw.data = subset.matrix.s2, project = 's2')

#-----------------------------------------------------------------------------
mito.genes <- grep(pattern = "^mt-", x = rownames(x = so.s1@data), value = TRUE)
percent.mito <- Matrix::colSums(so.s1@raw.data[mito.genes, ])/Matrix::colSums(so.s1@raw.data)
so.s1 <- AddMetaData(object = so.s1, metadata = percent.mito, col.name = "percent.mito")
#-----------------------------------------------------------------------------
mito.genes <- grep(pattern = "^mt-", x = rownames(x = so.s2@data), value = TRUE)
percent.mito <- Matrix::colSums(so.s2@raw.data[mito.genes, ])/Matrix::colSums(so.s2@raw.data)
so.s2 <- AddMetaData(object = so.s2, metadata = percent.mito, col.name = "percent.mito")
#-----------------------------------------------------------------------------

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

so.s2 <- ScaleData(object = so.s2, vars.to.regress = c("nUMI", "percent.mito"))
#-----------------------------------------------------------------------------
so.org <- so

name.s1 <- 'ctr'
name.s2 <- 'cko'

wd.org <- getwd()
system(paste('mkdir -p', name.s1))
setwd(name.s1)
so <- so.s1
save(so, file='so.rdb')
setwd(wd.org)

system(paste('mkdir -p', name.s2))
setwd(name.s2)
so <- so.s2
save(so, file='so.rdb')
setwd(wd.org)

#-----------------------------------------------------------------------------







