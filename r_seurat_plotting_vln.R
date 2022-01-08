#!/usr/bin/Rscript
#------------------------------------------------------------------------------------
# for Soo

args <- commandArgs(TRUE)
if(length(args) < 1) {
	cat("\nUsage: r_seurat_plotting_vln.R [SO_PREFIX]
	./r_seurat_plotting_vln.R clu_10\n\n")
	q()
}

so.prefix <- args[1]

file.genes <- '/data/abby/soo/soo_20180423_scrnaseq_jmjd3_mn_brc_thr/rep1_jmjd3_e12_mn_brc_thr_control/52_known_markers.gs.all'

#file.genes <- '/home/abby/gene.gs'

library(Seurat)


load('so.rdb')

#ident.ordered <- so@ident
#levels(ident.ordered) = c("2", "8", "10", "5", "1", "4", "12", "0", "3", "9", "11", "6", "7")
#so@ident <- ident.ordered

genes.features.plot <- sort(unique(read.table(file.genes, stringsAsFactors=F)$V1))

genes.in.the.list <- sort(rownames(so@raw.data)[which(rownames(so@raw.data) %in% genes.features.plot)])

n.pages <- ceiling(length(genes.in.the.list) / 9)

#so.prefix <- 'comp_rep2_control_vs_rep2_ko_clu10'
#-----------------------------------
dir.out <- paste(so.prefix,'_plot_vln_markers',sep="")
system(paste("rm -rf ", dir.out, "; mkdir ",dir.out,sep=""))

i<-1
	i.start <- 9 * ( i - 1 ) + 1
	i.end <- 9 * i
	part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
	part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
	png(file.png)
	VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
	dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
	png(file.png)
	VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
	dev.off()

i<-2
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-3
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-4
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-5
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-6
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-7
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-8
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-9
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()

q()
i<-10
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-11
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
i<-12
        i.start <- 9 * ( i - 1 ) + 1
        i.end <- 9 * i
        part.genes.in.the.list <- genes.in.the.list[c(i.start:i.end)]
        part.genes.in.the.list <- part.genes.in.the.list[!is.na(part.genes.in.the.list)]
        file.png <- paste(dir.out,"/wo_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
        dev.off()
        file.png <- paste(dir.out,"/wt_vln_markers_",i,".png",sep="")
        png(file.png)
        VlnPlot(object = so, point.size.use = 1, features.plot = part.genes.in.the.list )
        dev.off()
