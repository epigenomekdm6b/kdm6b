#!/usr/bin/Rscript

PercentAbove <- function(x, threshold){
	  return(length(x = x[x > threshold]) / length(x = x))
}

CountAbove <- function(x, threshold){
	  return(length(x = x[x > threshold]))
}


label.y2 <- c(
"Mnx1_co", "Mnx1_ko",
"Isl1_co", "Isl1_ko",
"Isl2_co", "Isl2_ko",
"Lhx3_co", "Lhx3_ko",
"Lhx4_co", "Lhx4_ko",
"Crabp1_co", "Crabp1_ko",
"Klf5_co", "Klf5_ko",
"Lix1_co", "Lix1_ko",
"Foxp1_co", "Foxp1_ko",
"Aldh1a2_co", "Aldh1a2_ko",
"Lhx1_co", "Lhx1_ko",
"Kcnip4_co", "Kcnip4_ko",
"Kitl_co", "Kitl_ko",
"Nos1_co", "Nos1_ko",
"Nkx2.2_co", "Nkx2.2_ko",
"Hoxc6_co", "Hoxc6_ko",
"Hoxc9_co", "Hoxc9_ko",
"Tubb3_co", "Tubb3_ko",
"Nefl_co", "Nefl_ko" 
)
label.y <- c(
"co_Mnx1", "ko_Mnx1",
"co_Isl1", "ko_Isl1",
"co_Isl2", "ko_Isl2",
"co_Lhx3", "ko_Lhx3",
"co_Lhx4", "ko_Lhx4",
"co_Crabp1", "ko_Crabp1",
"co_Klf5", "ko_Klf5",
"co_Lix1", "ko_Lix1",
"co_Foxp1", "ko_Foxp1",
"co_Aldh1a2", "ko_Aldh1a2",
"co_Lhx1", "ko_Lhx1",
"co_Kcnip4", "ko_Kcnip4",
"co_Kitl", "ko_Kitl",
"co_Nos1", "ko_Nos1",
"co_Nkx2.2", "ko_Nkx2.2",
"co_Hoxc6", "ko_Hoxc6",
"co_Hoxc9", "ko_Hoxc9",
"co_Tubb3", "ko_Tubb3",
"co_Nefl", "ko_Nefl" 
)
levels.gene.ordered <- c(
"Mnx1_CTRL", "Mnx1_KO",
"Isl1_CTRL", "Isl1_KO",
"Isl2_CTRL", "Isl2_KO",
"Lhx3_CTRL", "Lhx3_KO",
"Lhx4_CTRL", "Lhx4_KO",
"Crabp1_CTRL", "Crabp1_KO",
"Klf5_CTRL", "Klf5_KO",
"Lix1_CTRL", "Lix1_KO",
"Foxp1_CTRL", "Foxp1_KO",
"Aldh1a2_CTRL", "Aldh1a2_KO",
"Lhx1_CTRL", "Lhx1_KO",
"Kcnip4_CTRL", "Kcnip4_KO",
"Kitl_CTRL", "Kitl_KO",
"Nos1_CTRL", "Nos1_KO",
"Nkx2.2_CTRL", "Nkx2.2_KO",
"Hoxc6_CTRL", "Hoxc6_KO",
"Hoxc9_CTRL", "Hoxc9_KO",
"Tubb3_CTRL", "Tubb3_KO",
"Nefl_CTRL", "Nefl_KO" 
)


#function (object, grouping.var, genes.plot, gene.groups, cols.use = c("blue", 
    #"red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    #group.by, plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
#{

load('so.rdb')
library(Seurat)
library(tidyr)
library(dplyr)

file.genes <- '/data/abby/soo/soo_20180423_scrnaseq_jmjd3_mn_brc_thr/compare_jmjd3_e12_mn_brc_thr_rep2_control_vs_rep2_ko/seurat_comp_rep2_control_vs_ko/genes_for_drawing.gs'
genes.plot <- unique(read.table(file.genes, stringsAsFactors=F)$V1)


object=so
grouping.var = 'stim' 
#genes.plot, 
#gene.groups, 
cols.use = c("blue", "red")
col.min = -2.5
col.max = 2.5
dot.min = 0
dot.scale = 6
    #group.by, 
plot.legend = FALSE
do.return = FALSE
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

    data.to.plot <- data.to.plot %>% unite(id, c("id","class"))
#---------------------------------------------------------------------

	data.to.plot <- data.to.plot %>% separate(id, into=c("clu","class"), sep="_") %>% 
			#filter(clu==2 | clu==5 | clu==4 | clu==1 | clu==12 | clu==3 | clu==0 | clu==9 | clu==11 | clu==6 | clu==7) %>%
			filter(clu==2 | clu==8 | clu==5 | clu==4 | clu==1 | clu==12 | clu==3 | clu==0) %>%
			#*****************************************************************
			unite(id, c("clu","class"))

    	data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
        	value = expression, -c(cell, id))
#---------------------------------------------------------------------
	data.to.plot <- data.to.plot %>% separate(id,into=c("clu","class"), sep="_") %>% 
		mutate(id.clu=clu) %>% 
		unite("id",c("genes.plot","class"),sep="_")
                 		#cell clu        id expression id.clu
		#1 s1_AAACCTGGTAGCCTCG   1 Mnx1_CTRL  0.5614973      1
		#2 s1_AAACCTGGTCTAGTCA   7 Mnx1_CTRL  0.0000000      7
#--------------------------------------------------------
	data.to.plot <- data.to.plot %>% group_by(clu, id) %>%
		summarise(avg.exp = ExpMean(x=expression), pct.exp = PercentAbove(x=expression, threshold=0))
		# A tibble: 6 x 4
		# Groups:   clu [1]
  		#clu   id           avg.exp pct.exp
  		#<chr> <chr>          <dbl>   <dbl>
		#1 0     Aldh1a2_CTRL 0.00985  0.0114
#--------------------------------------------------------
	data.to.plot <- data.to.plot %>% mutate(avg.exp=scale(x=avg.exp)) %>% 
		mutate(avg.exp.scale=as.numeric(x=cut(x=MinMax(data=avg.exp, max=col.max, min=col.min), breaks=100)))
		# A tibble: 6 x 5
		# Groups:   clu [1]
  		#clu   id           avg.exp pct.exp avg.exp.scale
  		#<chr> <chr>          <dbl>   <dbl>         <dbl>
		#1 0     Aldh1a2_CTRL -0.638   0.0114             1
		#2 0     Aldh1a2_KO   -0.631   0.0202             1

    #data.to.plot <- data.to.plot %>% separate(col = id, into = c("ident1", 
        #"ident2"), sep = "_") %>% rowwise() %>% mutate(palette.use = "black", 
        #ptcolor = colorRampPalette(colors = c("grey", "black"))(20)[avg.exp.scale]) %>% 
        #unite("id", c("ident1", "ident2"), sep = "_")

#-------------------------------------------------
fun.ptcolor <- function(ident2,avg.exp.scale) {
#-------------------------------------------------
	ptcolor <- ifelse(ident2=="CTRL",
			colorRampPalette(colors = c("lightblue", "blue"))(100)[avg.exp.scale],
			colorRampPalette(colors = c("pink", "red"))(100)[avg.exp.scale]
	)
	ptcolor
}
#-------------------------------------------------

	data.to.plot <- data.to.plot %>% separate(col=id,into=c("ident1","ident2"),sep="_") %>% rowwise() %>% 
		mutate(palette.use=ifelse(ident2=="CTRL","blue","red"),
		       ptcolor=fun.ptcolor(ident2, avg.exp.scale)) %>%
        	unite("id", c("ident1", "ident2"), sep = "_")

	#mutate(gradebook, Pass.Fail = ifelse(grade > 60, "Pass", "Fail"))


	ordered <- rev(levels.gene.ordered)
	data.to.plot$id <- factor(data.to.plot$id, levels=ordered)
	data.to.plot$clu <- factor(x=data.to.plot$clu, levels=c(2, 8, 5, 4, 1, 12, 3, 0)) # for MNs

			#*****************************************************************
	scale.y.discrete <- rev(label.y2)

#2, 5, 4, 1, 12, 3, 0, 9, 11, 6, 7

    data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
    palette.use <- unique(x = data.to.plot$palette.use)




p <- ggplot(data = data.to.plot, mapping = aes(x = clu, y = id)) + 
	geom_point(mapping = aes(size = pct.exp, color = ptcolor)) +
	scale_radius(range = c(0, dot.scale)) + scale_color_identity() + 
	#scale_x_discrete(labels=c("C2","C5","C4","C1","C12","C3","C0","C9","C11","C6","C7"), position="top") +
	scale_x_discrete(labels=c( "MMC", "MMC.Sp8+", "HMC", "LMCm", "LMCl", "LMCd", "PGC.Isl1+", "PGC.Isl1-"), position="top") +
	scale_y_discrete(labels=scale.y.discrete) +
	theme_bw() +
	theme(axis.title=element_blank())

#pdf("out1.pdf"); p; dev.off()
png("out1.png"); p; dev.off()


