#!/usr/bin/Rscript
#------------------------------------------------------------------------------------
# for Soo
rm(list=ls())
load('so.rdb')
library(Seurat)
library(ggplot2)

	clu.name <- c(
		      "mmc",
		      "mmc.sp8.pos",
		      "hmc.or.mmc",
		      "hmc",
		      "lmcl.i.neg",
		      "lmcm.i.pos",
		      "lmc.digit",
		      "pgc.i.neg",
		      "pgc.i.pos",
		      "mn.nb",
		      "mn.p",
		      "in.v2",
		      "in.v3"
		      )

	clu.no <- c("2", "8", "10", "5", "1", "4", "12", "0", "3", "9", "11", "6", "7")

        ident <- levels(so@ident)
        features.plot <- c("Zeb2")
        data <- data.frame(FetchData(object=so,vars.all=features.plot))

ident.ordered <- so@ident
levels(ident.ordered) = c("2", "8", "10", "5", "1", "4", "12", "0", "3", "9", "11", "6", "7")
so@ident <- ident.ordered

        data$ident <- so@ident

        colnames(data) <- c("feature","ident")
	y.max <- ceiling(max(data$feature))

        p <- ggplot(data = data, mapping = aes(x = factor(x = ident), y = feature))
        p <- p + geom_violin(scale = "width",trim = FALSE, mapping = aes(fill = factor(x = ident)))
        #p <- p + scale_x_discrete(labels=c("xx2", "8", "10", "5", "1", "4", "12", "0", "3", "9", "11", "6", "7"))
        p <- p + scale_x_discrete(labels=clu.name)
        p <- p + guides(fill = guide_legend(title = NULL)) + xlab("Identity")
        p <- p + theme(legend.position = "none")
        p <- p + ggpubr::rotate_x_text()
        p <- p + ylim(0,y.max)
        p <- p + ggtitle(features.plot)

	png('tmp.png')
        p
	dev.off()
#------------------------------------------------


load('so.rdb')


	png(file.png)
	VlnPlot(object = so, point.size.use = 0, features.plot = part.genes.in.the.list )
	dev.off()
