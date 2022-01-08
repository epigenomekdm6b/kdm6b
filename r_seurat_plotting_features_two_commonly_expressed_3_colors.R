#!/home/groups/precepts/chhy/anaconda3/envs/jmjd3/bin/Rscript

#------------------------------------------------------------------
fun.ix.cell.common <- function(so,genes.target) {

        col.gene <- FetchData(object=so,vars.all=genes.target)

        cutoff <- 0.1
        #cutoff <- 0.05
        ix.gt.in.gene1 <- which(col.gene[,1] > cutoff)
        ix.gt.in.gene2 <- which(col.gene[,2] > cutoff)

        ix.cell.common <- intersect(ix.gt.in.gene1,ix.gt.in.gene2)
        return(ix.cell.common)
}

color.pattern.map <- function(x,col.pattern) {
        col <- colorRampPalette(colors=c("grey",col.pattern))(20)[x]
        return(col)
}

alpha.level <- function(x) {
        if(x==1) { alpha <- 0.05 }
        else { alpha <- 0.6 }
        return(alpha)
}

library(Seurat)
load("so.rdb")

#--------------------------------------------------------------------------
dir.out <- paste("fig3_merged_feature",sep="")
#system(paste("rm -rf ", dir.out, "; mkdir ",dir.out,sep=""))

#--------------------------------------------------------------------------
# abby
genes.target <- c("Slc18a3","Slc17a6")
#genes.target <- c("Isl1","Lhx3")
#genes.target <- c("Foxp1","Aldh1a2")
#--------------------------------------------------------------------------

colors.pattern <- c("red","dark green")

col.tsne <- GetCellEmbeddings(object=so,reduction.type="tsne")


#------------------------------------------------------------------
# for 1st gene
#------------------------------------------------------------------
	gene.target <- genes.target[1]
	color.pattern <- colors.pattern[1]

	col.gene <- FetchData(object=so,vars.all=gene.target)

	avg.exp.scale <- as.numeric(cut(col.gene[,1],breaks=20))

	ptcolor <- unlist(lapply(avg.exp.scale,color.pattern.map,col.pattern=color.pattern))

	alpha <- unlist(lapply(avg.exp.scale,alpha.level))

	data.1 <- data.frame(col.tsne, col.gene, avg.exp.scale, ptcolor, alpha)

#------------------------------------------------------------------
# for 2nd gene
#------------------------------------------------------------------
        gene.target <- genes.target[2]
        color.pattern <- colors.pattern[2]

        col.gene <- FetchData(object=so,vars.all=gene.target)

        avg.exp.scale <- as.numeric(cut(col.gene[,1],breaks=20))

        ptcolor <- unlist(lapply(avg.exp.scale,color.pattern.map,col.pattern=color.pattern))

        alpha <- unlist(lapply(avg.exp.scale,alpha.level))

        data.2 <- data.frame(col.tsne, col.gene, avg.exp.scale, ptcolor, alpha)


#------------------------------------------------------------------
data <- data.1
sc.1 <- data.1$avg.exp.scale
sc.2 <- data.2$avg.exp.scale

t <- rep(0,dim(data)[1])
df <- data.frame(t)
ix.1 <- which(sc.1>=sc.2)
ix.2 <- which(sc.1<sc.2)

df$col.1 <- data.1$ptcolor
df$col.2 <- data.2$ptcolor

df$ap.1 <- data.1$alpha
df$ap.2 <- data.2$alpha

df$sc.1 <- sc.1
df$sc.2 <- sc.2

#df[ix.1,"t"] <- "#FF99FF"
#df[ix.2,"t"] <- "#3399FF"
#df[ix.cell.common,"t"] <- "#FF0000"

rownames(df) <- rownames(data)
df$tSNE_1 <- data$tSNE_1
df$tSNE_2 <- data$tSNE_2

ptcolor <- c()
alpha <- c()
for(i in 1:dim(df)[1]) {
	sc.1 <- df[i,'sc.1']
	sc.2 <- df[i,'sc.2']

	col.1 <- as.character(df[i,'col.1'])
	col.2 <- as.character(df[i,'col.2'])

	ap.1 <- df[i,'ap.1']
	ap.2 <- df[i,'ap.2']

	ptcolor <- c(ptcolor, ifelse(sc.1>=sc.2, col.1, col.2))
	alpha <- c(alpha, ifelse(sc.1>=sc.2, ap.1, ap.2))
}
df$ptcolor <- as.character(ptcolor)
df$alpha <- alpha

ix.cell.common <- fun.ix.cell.common(so, genes.target)

#df[ix.cell.common,"ptcolor"] = "#FFFF00" ### Yellow
#df[ix.cell.common,"ptcolor"] = "#800080" ### Purple
df[ix.cell.common,"ptcolor"] = "#00FFFF" ### Cyan
#df[ix.cell.common,"ptcolor"] = "#006400" ### Dark green
df[ix.cell.common,"alpha"] = 0.6

df$ptcolor <- as.factor(df$ptcolor)

#------------------------------------------------------------------
label <- paste0(genes.target[1],"_",genes.target[2])
	
	p <- ggplot(df,aes(x=tSNE_1,y=tSNE_2)) +
		#geom_point(aes(colour=ptcolor, alpha=alpha), show.legend=F, size=0.5) +
		geom_point(aes(colour=ptcolor, alpha=alpha), show.legend=F, size=1) +
                theme_bw() +
                theme(
                        panel.grid = element_blank(),
                        #text = element_text(size=30),
                        plot.title = element_text(hjust = 0.5)
                ) +
		scale_color_identity() +
		ggtitle(label)


# png output
out.img <- paste0(dir.out,"/", "merged_",genes.target[1],"_",genes.target[2],".png")
png(out.img, bg='transparent'); print(p); dev.off()

# pdf output
out.img <- paste0(dir.out,"/", "merged_",genes.target[1],"_",genes.target[2],".pdf")
pdf(out.img, bg='transparent',width=3,height=3); print(p); dev.off()

#--------------------------------------------------------------------------


