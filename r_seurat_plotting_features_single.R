#!/home/groups/precepts/chhy/anaconda3/envs/jmjd3/bin/Rscript

args <- commandArgs(TRUE)
if(length(args) < 1) {
        cat("\nUsage: r_seurat_plotting_features_single.R [GS_FILE]
        ./r_seurat_plotting_features_single.R selected_markers.gs\n\n")
        q()
}

gs.file <- args[1]

so.prefix <- gsub(".gs","",basename(gs.file))

genes.target <- read.table(gs.file, stringsAsFactors=F)$V1

color.pattern.map <- function(x,col.pattern) {
        #col <- colorRampPalette(colors=c("grey",col.pattern))(20)[x]
        col <- colorRampPalette(colors=c("grey",col.pattern))(20)[x]
        return(col)
}

alpha.level <- function(x) {
        if(x==1) { alpha <- 0.3 }
        else { alpha <- 0.4 }
        return(alpha)
}

library(Seurat)
load("so.rdb")

dir.out <- paste0(so.prefix,"_feature")
system(paste("rm -rf ", dir.out, "; mkdir ",dir.out,sep=""))

#genes.target <- c("Cited")
#genes.target <- c("Npy","Ghrh","Sox14")
#colors.pattern <- c("blue","red","green")
color.pattern <- "blue"

col.tsne <- GetCellEmbeddings(object=so,reduction.type="tsne")

n.genes <- length(genes.target) 

for(i in 1:n.genes) {

	gene.target <- genes.target[i]

	col.gene <- FetchData(object=so,vars.all=gene.target)

	avg.exp.scale <- as.numeric(cut(col.gene[,1],breaks=20))

	ptcolor <- unlist(lapply(avg.exp.scale,color.pattern.map,col.pattern=color.pattern))

	alpha <- unlist(lapply(avg.exp.scale,alpha.level))

	data <- data.frame(col.tsne, col.gene, avg.exp.scale, ptcolor, alpha)

	p <- ggplot(data,aes(x=tSNE_1,y=tSNE_2)) +
		#geom_point(aes(colour=ptcolor, alpha=alpha), show.legend=F, size=0.5) +
		geom_point(aes(colour=ptcolor, alpha=alpha), show.legend=F, size=2) +
		#geom_point(aes(colour=ptcolor, alpha=alpha), show.legend=F, size=1) +
		scale_color_identity() +
		ggtitle(gene.target)
	#p <- p + theme_bw()
	p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.line = element_line(colour = "black"))

	# png output
	img.out <- paste0(dir.out,'/plot_features_',gene.target,'.png'); 
	png(img.out,width=300,height=300,bg='transparent'); print(p); dev.off()

	# pdf output
	img.out <- paste0(dir.out,'/plot_features_',gene.target,'.pdf'); 
	pdf(img.out,width=3,height=3,bg='transparent'); print(p); dev.off()

}




