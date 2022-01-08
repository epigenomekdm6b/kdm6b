#!/home/groups/precepts/chhy/anaconda3/envs/jmjd3/bin/Rscript

args <- commandArgs(TRUE)
if(length(args) < 2) {
	cat("\nUsage: r_seurat_plotting_vln_for_manuscript_stacking.R [GS_FILE] [no|new]
	./r_seurat_plotting_vln_for_manuscript_stacking.R selected_markers.gs no|new\n\n")
	q()
}

file.genes <- args[1]
new.label <- args[2]

load('so.rdb')

so.prefix <- gsub(".gs","",basename(file.genes))

library(Seurat)

fun.clu.name <- function(clu.profile, clu.target.id) {
	clu.target.id <- as.character(clu.target.id)
	hs <- cluster$name
	names(hs) <- cluster$id
	return(as.character(unname(hs[clu.target.id])))
}

if( new.label == "new" ) {
	#***************************************************************************
	# checking cluster_file
	cluster.file <- 'cluster.profile'
	if(!file.exists(cluster.file)) { print("cluster.prorile file required"); q(); }
	cluster <- read.table(cluster.file, stringsAsFactors=T, row.names=1, header=T, sep="\t")

	ident.named <- unlist(lapply(unname(so@ident),fun.clu.name,clu.profile=cluster))
	ident.named <- factor(ident.named, level=cluster$name[cluster$order])

	#------------------------------------------------------
	# New labeling....
	so@ident <- ident.named 
	names(so@ident) <- colnames(so@data)

	cols.use <- paste0("#",cluster$color)
}

#***************************************************************************
genes.in.the.list <- read.table(file.genes, stringsAsFactors=F)$V1

#------------------------------------------------------------------------------------
dir.out <- paste(so.prefix,"_vln",sep="")
system(paste("rm -rf ", dir.out, "; mkdir ",dir.out,sep=""))

#--------------------------------------
draw <- function(gene.target) {
#--------------------------------------
        p <- VlnPlot(object = so, point.size.use = NA, features.plot = c(gene.target), 
		size.x.use = 0, 
		size.y.use = 1, 
		size.title.use = 7,
		x.lab.rot=T, cols.use=cols.use
	) + theme(plot.margin = unit(c(0.2,0.2,0,0), "cm")) + theme(plot.title=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_text(size=8,family="Courier"))

        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

	return(p)
}

#--------------------------------------
# returned ggplot object2
#--------------------------------------
plotlist <- lapply(genes.in.the.list, draw)
height <- length(genes.in.the.list)

#--------------------------------------
# print individual plot
#--------------------------------------
	for( i in 1:length(genes.in.the.list) ) {
		p <- plotlist[[i]]
		gene <- genes.in.the.list[[i]]

        	file.img <- paste(dir.out,"/","individual_vln_",gene,".pdf",sep="")
        	pdf(file=file.img, bg='transparent', width=3, height=1); print(p); dev.off()

        	file.img <- paste(dir.out,"/","individual_vln_",gene,".png",sep="")
        	png(file=file.img, bg='transparent', width=300, height=100); print(p); dev.off()
	}


#--------------------------------------
# print stacked_plot
#--------------------------------------
	p <- plot_grid(plotlist=plotlist,nrow=length(plotlist),ncol=1,align="v")

	# pdf output
        file.img <- paste(dir.out,"/","vln_stacked_",so.prefix,".pdf",sep="")
        pdf(file=file.img, bg='transparent', width=3, height=height)
	print(p)
        dev.off()

	# png output
        file.img <- paste(dir.out,"/","vln_stacked_",so.prefix,".png",sep="")
	png(file.img, bg='transparent', width=300, height=100*height)
	print(p)
        dev.off()
#--------------------------------------
q()
