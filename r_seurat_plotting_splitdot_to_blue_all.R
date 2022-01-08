#!/usr/bin/Rscript

PercentAbove <- function(x, threshold){ return(length(x = x[x > threshold]) / length(x = x)) }
CountAbove <- function(x, threshold){ return(length(x = x[x > threshold])) }

load('so.rdb')
library(Seurat)
library(tidyr)
library(dplyr)

#---------------------------------------------------------------------
#file.genes <- '/data/abby/soo/soo_20180423_scrnaseq_jmjd3_mn_brc_thr/compare_jmjd3_e12_mn_brc_thr_rep2_control_vs_rep2_ko/seurat_comp_rep2_control_vs_ko/genes_for_drawing.gs'
#file.genes <- 'marker_dotplot_mn_26.gs'
#file.genes <- 'marker_dotplot_etc_17.gs'
#file.genes <- 'marker_dotplot_all_43.gs'

file.genes <- 'marker_dotplot_mn_new.gs'
#file.genes <- 'marker_dotplot_etc_new.gs'

#pMN	11
#nbMN	9
#V2IN	6
#V3IN	7
#ND	10

levels.clu.class.ordered <- c(
			      "11_CTRL", "11_KO",
			      "9_CTRL", "9_KO",
                              "2_CTRL", "2_KO",
                              "8_CTRL", "8_KO",
                              "5_CTRL", "5_KO",
                              "4_CTRL", "4_KO",
                              "1_CTRL", "1_KO",
                              "12_CTRL", "12_KO",
                              "3_CTRL", "3_KO",
                              "0_CTRL", "0_KO",
			      "6_CTRL", "6_KO",
			      "7_CTRL", "7_KO",
			      "10_CTRL", "10_KO"
			      )

label.clu.class.ordered <- c(
			     "pMN.ctr", "pMN.cko",
			     "nbMN.ctr", "nbMN.cko",
                             "MMC.ctr", "MMC.cko",
                             "MMC.Sp8+.ctr", "MMC.Sp8+.cko",
                             "HMC.ctr", "HMC.cko",
                             "LMCm.ctr", "LMCm.cko",
                             "LMCl.ctr", "LMCl.cko",
                             "LMCd.ctr", "LMCd.cko",
                             "PGC.Isl1+.ctr", "PGC.Isl1+.cko",
                             "PGC.Isl1-.ctr", "PGC.Isl1-.cko",
			     "V2IN.ctr", "V2IN.cko",
			     "V3IN.ctr", "V3IN.cko",
			     "ND.ctr", "ND.cko"
			     )
#---------------------------------------------------------------------


genes.plot <- unique(read.table(file.genes, stringsAsFactors=F)$V1)
levels.gene.ordered <- genes.plot

object=so
col.min = -2.5
col.max = 2.5
dot.min = 0
dot.scale = 6
    #group.by, 
x.lab.rot = FALSE

#---------------------------------------------------------------------
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
    data.to.plot$cell <- rownames(x = data.to.plot)
    data.to.plot$id <- object@ident

    ix.s1 <- grep('s1_',data.to.plot$cell)
    ix.s2 <- grep('s2_',data.to.plot$cell)
    class <- vector("list", length=length(data.to.plot$id))
    class[ix.s1] <- "CTRL"
    class[ix.s2] <- "KO"
    data.to.plot$class <- unlist(class)

    #--------------------------------------------------------------------------
    # To filter in target clusters: 
    #					2,8,5,4,1,12,3,0 ordered for MNs
    #--------------------------------------------------------------------------
	#data.to.plot <- data.to.plot %>%
			#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			#filter(id==2 | id==8 | id==5 | id==4 | id==1 | id==12 | id==3 | id==0) # for MNs
			#filter(id==11 | id==9 | id==6 | id==7 | id==10 ) # for ETC
			#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #--------------------------------------------------------------------------

    #data.to.plot <- data.to.plot %>% unite(id, c("id","class"))
    tmp <- data.to.plot %>% unite(id, c("id","class"))
    tmp <- data.to.plot %>% gather(key=)
    	 tmp <- data.to.plot %>% gather(key = genes.plot, value = expression, -c(cell, id, class)) %>% 
		 	unite(clu.class, c("id","class"),sep="_")

	tmp1 <- tmp %>% group_by(clu.class,genes.plot) %>%
	summarise(avg.exp = ExpMean(x=expression), pct.exp = PercentAbove(x=expression, threshold=0))

        tmp1 <- tmp1 %>% mutate(avg.exp=scale(x=avg.exp)) %>%
                mutate(avg.exp.scale=as.numeric(x=cut(x=MinMax(data=avg.exp, max=col.max, min=col.min), breaks=100)))

	#-------------------------------------------------
	# for only Nkx2.2 
	tmp1$genes.plot <- gsub("\\.","-",tmp1$genes.plot)
	#-------------------------------------------------


#-------------------------------------------------
fun.ptcolor <- function(ident2,avg.exp.scale) {
#-------------------------------------------------
        ptcolor <- ifelse(ident2=="CTRL",
                        colorRampPalette(colors = c("lightblue", "blue"))(100)[avg.exp.scale],
                        #colorRampPalette(colors = c("yellow", "green"))(100)[avg.exp.scale],
                        colorRampPalette(colors = c("pink", "red"))(100)[avg.exp.scale]
        )
        ptcolor
}
#-------------------------------------------------

        tmp1 <- tmp1 %>% separate(col=clu.class,into=c("ident1","ident2"),sep="_") %>% rowwise() %>%
                mutate(palette.use=ifelse(ident2=="CTRL","blue","red"),
                       ptcolor=fun.ptcolor(ident2, avg.exp.scale)) %>%
                unite("clu.class", c("ident1", "ident2"), sep = "_")


	ordered <- rev(levels.gene.ordered)
	tmp1$genes.plot <- factor(tmp1$genes.plot, levels=ordered)
	tmp1$clu.class <- factor(x=tmp1$clu.class, levels=levels.clu.class.ordered) # for MNs

#---------------------------------------------------------------------
# temporary
p <- ggplot(data = tmp1, mapping = aes(x = clu.class, y = genes.plot)) + 
	geom_point(mapping = aes(size = pct.exp, color = ptcolor)) +
	scale_radius(range = c(0, dot.scale)) + scale_color_identity() + 
	scale_x_discrete(labels=label.clu.class.ordered, position="top") +
	#scale_y_discrete(labels=scale.y.discrete) +
	theme_bw() +
	theme(axis.title=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0, hjust=0))


pdf("plot_cko_dotplot_mn_new.pdf"); p; dev.off()
#pdf("plot_cko_dotplot_etc_new.pdf"); p; dev.off()
#pdf("plot_cko_dotplot_all.pdf"); p; dev.off()
#png("out1_horiz_etc_1.png"); p; dev.off()
#---------------------------------------------------------------------
