# Final figures into panels for publication

# Library stuff
library(ggplot2)
library(reshape)
library(plyr)
library(gtools)
library(gdata)
library(grDevices)
library(ggbiplot)
library(Rmisc)
library(grid)
library(cowplot)
library(ggfortify)
library(infotheo)

load("~/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/redefine_piles/.RData")

	reads2rpkm <- function (m_reads, tr_length, total_reads) {
	  
	  RPK <- (m_reads / (tr_length / 1000)) # mapped reads / transcript length in kb
	  RPKM <- (RPK / (total_reads /  1000000)) # RPK / total num reads in million
	  
	  return (RPKM)
	}


df.gene_counts.DNA.passages <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/data/Jan15/gene_counts_passages/axam_Ht0_Hp0_Hp10_Hp21.bed", header=TRUE)
df.gene_counts.DNA.passages$Chromsome = gsub(df.gene_counts.DNA.passages$Chromsome, pattern="\\|\\S+", replacement="")
colnames(df.gene_counts.DNA.passages) <- c("Chromosome","Start","End","Geneid","Score","Strand","Ax_ama" ,"Sp_ama","p2","p10","p21")
gpkm.passages <- df.gene_counts.DNA.passages
gpkm.passages$Chromosome <- factor(gpkm.passages$Chromosome, levels = mixedsort(unique(gpkm.passages$Chromosome)), ordered = TRUE)

T <- sum(gpkm.passages$Ax_ama)
gpkm.passages$Ax_ama<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[7]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$Ax_ama))})
T <- sum(gpkm.passages$Sp_ama)
gpkm.passages$Sp_ama<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[8]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$Sp_ama))})
T <- sum(gpkm.passages$p2)
gpkm.passages$p2<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[9]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p2))})
T <- sum(gpkm.passages$p10)
gpkm.passages$p10<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[10]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p10))})
T <- sum(gpkm.passages$p21)
gpkm.passages$p21<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[11]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p21))})
T <- sum(gpkm.passages$p20)
gpkm.passages$p20<-gpkm.passages$p21

gg_color_hue <- function(n) { 
	hues = seq(15, 375, length = n + 1) 
	hcl(h = hues, l = 65, c = 100)[1:n] 
}

gpkm.passages.m <- melt(gpkm.passages, measure.vars = c("Ax_ama", "Sp_ama", "p2", "p10", "p21", "p20"))

#################### Remeber to put back into supplementary the full boxplot for all chromosomes in all the passags



	gpkm = readRDS("~/projects/kinetoplastids/Ldonovani/data/R/gpkm.RDS")
	rpkm = readRDS("~/projects/kinetoplastids/Ldonovani/data/R/rpkm.RDS")
	gpkm.clones <- gpkm
	colnames(gpkm.clones) <- c("Gene","Chromosome","Start","End","Strand","Ax_ama","Sp_ama","CL1" ,"CL3","CL4","CL6","CL7", "CL8", "CL9", "CL10", "OSp_ama")
	gpkm.clones$Chromosome = gsub(gpkm.clones$Chromosome, pattern="(.*)\\|\\S+", replacement="")
	gpkm.clones$Chromosome <- factor(gpkm.clones$Chromosome, levels = mixedsort(unique(gpkm.clones$Chromosome)), ordered = TRUE)

	# Match clones and passage genes for plotting p20 together
	passage.genes.intersect <- subset(gpkm.passages, Geneid %in% gpkm.clones$Gene)
	gpkm.clones$p21 <- passage.genes.intersect[order(array(passage.genes.intersect$Geneid)), "p21"]
	gpkm.clones.m <- melt(gpkm.clones, measure.vars = c("Ax_ama", "Sp_ama", "OSp_ama", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10", "p21"))
	gpkm.clones.m$variable <- factor(gpkm.clones.m$variable, levels = c("Ax_ama", "Sp_ama", "OSp_ama", "p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10"), ordered = TRUE)

	# Get fucking medians per chromsome?
	gpkm.clones.m.medians <- ddply(gpkm.clones.m, .(variable, Chromosome), function(df) {
			(median(df$value) / median(subset(gpkm.clones.m, Chromosome == "34" & variable == unique(df$variable))$value))*2
		} )

	gpkm.clones.m2<- merge(gpkm.clones.m, gpkm.clones.m.medians, by = c("variable", "Chromosome"))


	df.sites = read.table("~/projects/kinetoplastids/Ldonovani/data/R/merged.wp20.grey.txt")
	colnames(df.sites) = c("chr", "position", "CL10", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "p20")
	df.sites$npos  = seq(from=1, to=dim(df.sites)[1])
	df.sites[which(df.sites$chr == "3"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9", "p20")] <- 5
	df.sites[which(df.sites$chr == "8"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9")] <- 5
	df.sites[which(df.sites$chr == "14"),c("CL1", "CL8")] <- 5
	df.sites[which(df.sites$chr == "15"),c("CL8")] <- 5
	df.sites[which(df.sites$chr == "16"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9")] <- 5
	df.sites[which(df.sites$chr == "20"),c("CL1", "CL8")] <- 5
	df.sites[which(df.sites$chr == "23"),c("CL1", "CL8")] <- 5
	df.sites[which(df.sites$chr == "32"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9", "p20")] <- 5
	df.sites.melt = melt(df.sites, id.vars=c("chr", "position", "npos"))
	df.sites.melt$chr <- factor(df.sites.melt$chr, levels = sort(unique(df.sites.melt$chr)), ordered = TRUE)

	gpkm.clones <- gpkm
	colnames(gpkm.clones) <- c("Gene","Chromosome","Start","End","Strand","Ax_ama","Sp_ama","CL1" ,"CL3","CL4","CL6","CL7", "CL8", "CL9", "CL10", "OSp_ama")
	gpkm.clones$Chromosome = gsub(gpkm.clones$Chromosome, pattern="(.*)\\|F.*", replacement="\\1")
	gpkm.clones$Chromosome <- factor(gpkm.clones$Chromosome, levels = mixedsort(unique(gpkm.clones$Chromosome)), ordered = TRUE)

    df.gene_counts.DNA.passages <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/data/Jan15/gene_counts_passages/axam_Ht0_Hp0_Hp10_Hp21.bed", header=TRUE)

	df.gene_counts.DNA.passages$Chromsome = gsub(df.gene_counts.DNA.passages$Chromsome, pattern="\\|\\S+", replacement="")
	colnames(df.gene_counts.DNA.passages) <- c("Chromosome","Start","End","Geneid","Score","Strand","Ax_ama" ,"Sp_ama","p2","p10","p21")
	gpkm.passages <- df.gene_counts.DNA.passages
	gpkm.passages$Chromosome <- factor(gpkm.passages$Chromosome, levels = mixedsort(unique(gpkm.passages$Chromosome)), ordered = TRUE)

	T <- sum(gpkm.passages$Ax_ama)
	gpkm.passages$Ax_ama<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[7]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$Ax_ama))})
	T <- sum(gpkm.passages$Sp_ama)
	gpkm.passages$Sp_ama<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[8]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$Sp_ama))})
	T <- sum(gpkm.passages$p2)
	gpkm.passages$p2<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[9]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p2))})
	T <- sum(gpkm.passages$p10)
	gpkm.passages$p10<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[10]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p10))})
	T <- sum(gpkm.passages$p21)
	gpkm.passages$p21<-apply(gpkm.passages,1, function(x){ reads2rpkm(as.numeric(x[11]),(as.numeric(x[3])-as.numeric(x[2])),sum(gpkm.passages$p21))})
	T <- sum(gpkm.passages$p20)
	gpkm.passages$p20<-gpkm.passages$p21

	gpkm.passages.m <- melt(gpkm.passages, measure.vars = c("Ax_ama", "Sp_ama", "p2", "p10", "p21", "p20"))

	# Match clones and passage genes for plotting p20 together
	passage.genes.intersect <- subset(gpkm.passages, Geneid %in% gpkm.clones$Gene)
	gpkm.clones$p21 <- passage.genes.intersect[order(array(passage.genes.intersect$Geneid)), "p21"]
	gpkm.clones.m <- melt(gpkm.clones, measure.vars = c("Ax_ama", "Sp_ama", "OSp_ama", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10", "p21"))
	gpkm.clones.m$variable <- factor(gpkm.clones.m$variable, levels = c("Ax_ama", "Sp_ama", "OSp_ama", "p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10"), ordered = TRUE)

	# Get fucking medians per chromsome?
	gpkm.clones.m.medians <- ddply(gpkm.clones.m, .(variable, Chromosome), function(df) {
			(median(df$value) / median(subset(gpkm.clones.m, Chromosome == "34" & variable == unique(df$variable))$value))*2
		} )
	
	#df.passages.piles <- read.table("~/projects/kinetoplastids/Ldonovani/data/R/Hp0_Hp10_Hp21_intersect_LDBPK.bed")
    df.passages.piles <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/analysis/SNP/new_passages_smiley/merged_new_piles.formatted4R.txt")

	colnames(df.passages.piles) <- c("chromosome", "start", "end", "sample","frequency")

	df.clones.piles <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/analysis/SNP_new/DNA/DNAseq_2014/pileup_smiley/DNA_clones.Rpile.txt")
	colnames(df.clones.piles) <- c("chromosome", "start", "end", "sample","frequency")
	colnames(df.clones.piles) <- c("chromosome", "start", "frequency", "sample","depth")

	df.t <- subset(df.passages.piles, !(chromosome %in% c("MC", "MI")) & sample == "Hp21")[,c("chromosome", "start", "sample", "frequency")]
	df.t$sample <- gsub(df.t$sample, pattern = "H", replacement = "")
	df.t2 <- subset(df.clones.piles, !(chromosome %in% c("MC", "MI")) & !(sample %in% c("DNA_0", "DNA_AX")))[,c("chromosome", "start", "sample", "frequency")]
	df.t2$sample <- gsub(df.t2$sample, pattern = "DNA_", replacement = "")
	df.pclones.piles <- rbind(df.t, df.t2)
	colnames(gpkm.clones.m.medians) <- c("sample", "chromosome", "somy")
	df.pclones.piles.medians <- merge(subset(gpkm.clones.m.medians, !(sample %in% c("Ax_ama", "Sp_ama", "OSp_ama")) & !(chromosome %in% c("MC", "MI"))), df.pclones.piles)
	df.pclones.somy = unique(df.pclones.piles.medians[,c("sample", "chromosome", "somy")])
	
	#
	df.pclones.tiles.count <- 	
	ddply(df.pclones.piles, .(chromosome, sample), function(d){ 
		print(unique(array(d$sample)))
		val = 0
		subpositions = unique(subset(df.sites, chr == unique(d$chromosome))$position)
		if(length(subpositions) > 0){
			val = length(unique(subset(d, start %in% subpositions)$start))[1]
		}
		return(val)
	})
	colnames(df.pclones.tiles.count) <- c("chromosome", "sample", "count")

	df.pclones.tiles.count.somy = merge(df.pclones.somy, df.pclones.tiles.count, by = c("chromosome", "sample"))

	df.pclones.tiles.count.somy$sample <- factor(df.pclones.tiles.count.somy$sample, levels = rev(c("p21","CL1", "CL8", "CL4", "CL3", "CL6", "CL7","CL9", "CL10")), ordered = TRUE)

	newvar = array(df.pclones.piles.medians$sample)
	newvarp = which(df.pclones.piles.medians$sample == "p21")
	newvar[newvarp] <- "p20"
	df.pclones.piles.medians$sample <- factor(newvar, levels = c("p20","CL1", "CL8", "CL4", "CL3", "CL6", "CL7","CL9", "CL10"), ordered = TRUE)
	df.sites.melt$variable <- factor(df.sites.melt$variable, levels = rev(c("p20","CL1", "CL8", "CL4", "CL3", "CL6", "CL7","CL9", "CL10")), ordered = TRUE)
	pcolors = colorRampPalette(c("blue", "green", "yellow", "red", "grey"))(5)
	i = 1
	plots <- list()


	labels=c("-","A","C","G","T")

		tmplot1 <- ggplot(subset(df.sites.melt, chr == "3"), aes(y=factor(variable), x=npos)) + 
		geom_tile(aes(fill=factor(value), height = .7)) + 
		scale_fill_manual(values=pcolors, labels=c("A","C","G","T", "none")) +
		scale_x_continuous(expand = c(0, 0)) + 
		facet_grid(.~chr, scales = "free") +
		#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
		labs(title = "", y = "", x = "", fill = "Dominant alelle") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(strip.background = element_blank(), plot.margin = unit(c(-0.3,0,-0.3,0), "cm"), 
			legend.direction = 'horizontal', legend.position = 'top') + 
			guides(fill=guide_legend(title = "Dominant allele", nrow=1,byrow=FALSE))

		tmplot2 <- ggplot(subset(df.pclones.piles.medians, chromosome == "3"), aes(frequency, fill = somy)) + 
		geom_density() + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		facet_wrap(~sample, ncol = 1 ) + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="", fill = "somy") + 
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(),
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank(),
			legend.direction = 'horizontal', legend.position = 'top') + 
			guides(fill=guide_legend(title = "somy",nrow=1,byrow=FALSE))


	
	tmplot3 <- ggplot(subset(df.pclones.tiles.count.somy, chromosome == "26"), aes(sample, count)) + 
	geom_bar(stat = "identity", width = .6, fill = 'white', colour = 'black') + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		coord_flip() + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="") + 
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(),
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank(),
			legend.direction = 'horizontal', legend.position = 'top') + 
			guides(fill=guide_legend(title = "somy",nrow=1,byrow=FALSE))


	i = 1
	#			 1    3    5    7    9     11    13    15    17    19    21    23			
	for (c in c("3", "5", "8", "9", "14", "15", "16", "20", "23", "26", "31", "32")) {  

		plots[[i]] <- ggplot(subset(df.sites.melt, chr == c), aes(y=variable, x=npos)) + 
		geom_tile(aes(fill=factor(value), height = .7)) + 
		scale_fill_manual(values=pcolors, labels=c("A","C","G","T"), breaks = c(0:3)) +
		scale_x_continuous(expand = c(0, 0)) + 
		facet_grid(.~chr, scales = "free") +
		#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
		labs(title = "", y = "", x = "", fill = "") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(legend.position = "none", strip.background = element_blank(), plot.margin = unit(c(-0.3,0,-0.3,0), "cm"))
		#panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(1,"lines")), legend.position = "none")
		#theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
		#science_theme2
		#plots[[i+1]] <- ggplot(subset(df.pclones.piles.medians, chromosome == c), aes(frequency, fill = somy)) + 
		#geom_density() + 
		#scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		#facet_wrap(~sample, ncol = 1 ) + 
		##labs(x = "", y="", fill = "Somy") + 
		#labs(x = "", y="", fill = "") + 
		#theme_bw(base_size = 8, base_family = "ArialMT") +
		#theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
		#	plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
		#	axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
		#	strip.text.x = element_blank(), legend.position = "none",
        #	axis.text.y=element_blank(),
        #	axis.ticks.y=element_blank())

		plots[[i+1]] <- ggplot(subset(df.pclones.tiles.count.somy, chromosome == c), aes(sample, count, fill = somy)) + 
		geom_bar(stat = "identity", width=.7) + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		coord_flip() + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="", fill = "somy") + 
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(), legend.position = "none",
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank())
		i = i + 2
	}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

  leg1 <- g_legend(tmplot1)
  leg2 <- g_legend(tmplot2)


xp1 = 0.25
yp1 = 0.31
yp2 = 0.325
wc1 = 0.18
diffp = 0.03
wc2 = 0.03
hc1 = 0.30
hc2 = 0.25
hc3 = 0.29


	test.p = read.table("~/projects/kinetoplastids/Ldonovani/analysis/SNP/new_passages_smiley/fixed_new_merged.txt")

	colnames(test.p) <- c("Chromosome", "Position", "Frequency", "Sample", "Depth")
	chrv <- array(test.p$Chromosome)
	test.p$Chromosome <- factor(chrv, levels=mixedsort(unique(chrv)), ordered = TRUE)
	test.p$Sample <- array(test.p$Sample)
	test.p[which(test.p$Sample == "Ht0"), "Sample"] <- "Splenic amastigote"
	test.p[which(test.p$Sample == "Hp0"), "Sample"] <- "p2"
	test.p[which(test.p$Sample == "Hp10"), "Sample"] <- "p10"
	test.p[which(test.p$Sample == "Hp21"), "Sample"] <- "p20"
	test.p$Sample <- factor(test.p$Sample, levels = c("Splenic amastigote", "p2", "p10", "p20"), ordered = TRUE)
	
	f4b <- ggplot(subset(test.p, Chromosome %in% c("5", "9", "14", "15", "20", "23", "26", "31")), aes(x = Frequency, colour = Sample)) + 
	geom_density() + 
	facet_wrap(~Chromosome, ncol=8) + 
	labs(x = "Frequency", y = "Density", title="") + 
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	science_theme2 +
	theme(
		legend.margin = unit(c(-0.6,0,0,0), "lines"),
    	axis.text.x=element_blank(), axis.ticks.x=element_blank(),  
	    legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') + 
		guides(fill=guide_legend(nrow=1,byrow=FALSE))


	postscript("~/plotsave/new_fig4_with_bars.eps", width = 8, height = 7, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")

	ggdraw() +
	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +
	draw_plot(plots[[3]],  0, 	yp1*2, wc1, hc1) +
	draw_plot(plots[[7]], xp1, 	yp1*2, wc1, hc1) +
	draw_plot(plots[[9]], xp1*2, yp1*2, wc1, hc1) +
	draw_plot(plots[[11]],xp1*3, yp1*2, wc1, hc1) +

	draw_plot(plots[[15]], 0, 	yp1, wc1, hc1) +
	draw_plot(plots[[17]], xp1,	yp1, wc1, hc1) +
	draw_plot(plots[[19]], xp1*2,	yp1, wc1, hc1) +
	draw_plot(plots[[21]], xp1*3,	yp1, wc1, hc1) +

	draw_plot(plots[[4]], xp1-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +	
	draw_plot(plots[[8]], (xp1*2)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +
	draw_plot(plots[[10]], (xp1*3)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +
	draw_plot(plots[[12]], (xp1*4)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +

	draw_plot(plots[[16]], xp1-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[18]], (xp1*2)-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[20]], (xp1*3)-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[22]], (xp1*4)-diffp, yp2+0.01, wc2, hc2-0.007) +

	draw_plot(f4b, 0, 0, 1, 0.26) +
  	draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.26), size = 16)

    dev.off()
    
    
xp1 = 0.25
yp1 = 0.31
yp2 = 0.32
wc1 = 0.18
diffp = 0.04
wc2 = 0.03
hc1 = 0.30
hc2 = 0.25
hc3 = 0.29


	postscript("~/plotsave/new_figs6new_with_bars.eps", width = 8, height = 7, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	
	ggdraw() +
	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +
	draw_plot(plots[[1]],		0,	yp1*2,	wc1,		hc1) +
	draw_plot(plots[[3]],	xp1*1,	yp1*2,	wc1,		hc1) +
	draw_plot(plots[[5]],	xp1*2, 	yp1*2,	wc1, 	hc1) +
	draw_plot(plots[[7]], 	xp1*3,	yp1*2,	wc1, 	hc1) +

	draw_plot(plots[[9]], 		0, 	yp1,	wc1, 	hc1) +
	draw_plot(plots[[11]], 	xp1*1,	yp1, 	wc1, 	hc1) +
	draw_plot(plots[[13]], 	xp1*2, 	yp1, 	wc1, 	hc1) +
	draw_plot(plots[[15]], 	xp1*3,	yp1, 	wc1, 	hc1) +

	draw_plot(plots[[17]], 		0, 		0,	wc1, 	hc1) +
	draw_plot(plots[[19]], 	xp1*1,		0,	wc1, 	hc1) +
	draw_plot(plots[[21]], 	xp1*2, 		0,	wc1, 	hc1) +
	draw_plot(plots[[23]], 	xp1*3,		0,	wc1, 	hc1) +

	draw_plot(plots[[2]], xp1-diffp, yp2*2, wc2, hc2) +	
	draw_plot(plots[[4]], (xp1*2)-diffp, yp2*2, wc2, hc2) +
	draw_plot(plots[[6]], (xp1*3)-diffp, yp2*2, wc2, hc2) +
	draw_plot(plots[[8]], (xp1*4)-diffp, yp2*2, wc2, hc2) +

	draw_plot(plots[[10]], xp1-diffp, yp2+0.01, wc2, hc2) +
	draw_plot(plots[[12]], (xp1*2)-diffp, yp2+0.01, wc2, hc2) +
	draw_plot(plots[[14]], (xp1*3)-diffp, yp2+0.01, wc2, hc2) +
	draw_plot(plots[[16]], (xp1*4)-diffp, yp2+0.01, wc2, hc2) +

	draw_plot(plots[[18]], xp1-diffp, 0.02, wc2, hc2) +	
	draw_plot(plots[[20]], (xp1*2)-diffp, 0.02, wc2, hc2) +
	draw_plot(plots[[22]], (xp1*3)-diffp, 0.02, wc2, hc2) +
	draw_plot(plots[[24]], (xp1*4)-diffp, 0.02, wc2, hc2)

dev.off()




stop("Trololo")



#### Another option





	i = 1
	#			 1    3    5    7    9     11    13    15    17    19    21    23			
	for (c in c("3", "5", "8", "9", "14", "15", "16", "20", "23", "26", "31", "32")) {  

		plots[[i]] <- ggplot(subset(df.sites.melt, chr == c), aes(y=variable, x=npos)) + 
		geom_tile(aes(fill=factor(value), height = .7)) + 
		scale_fill_manual(values=pcolors, labels=c("A","C","G","T"), breaks = c(0:3)) +
		scale_x_continuous(expand = c(0, 0)) + 
		facet_grid(.~chr, scales = "free") +
		#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
		labs(title = "", y = "", x = "", fill = "") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(legend.position = "none", strip.background = element_blank(), plot.margin = unit(c(-0.3,0,-0.3,0), "cm"))
		#panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(1,"lines")), legend.position = "none")
		#theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
		#science_theme2
		#plots[[i+1]] <- ggplot(subset(df.pclones.piles.medians, chromosome == c), aes(frequency, fill = somy)) + 
		#geom_density() + 
		#scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		#facet_wrap(~sample, ncol = 1 ) + 
		##labs(x = "", y="", fill = "Somy") + 
		#labs(x = "", y="", fill = "") + 
		#theme_bw(base_size = 8, base_family = "ArialMT") +
		#theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
		#	plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
		#	axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
		#	strip.text.x = element_blank(), legend.position = "none",
        #	axis.text.y=element_blank(),
        #	axis.ticks.y=element_blank())

		plots[[i+1]] <- ggplot(subset(df.pclones.tiles.count.somy, chromosome == c), aes(sample, count, fill = somy)) + 
		geom_bar(stat = "identity", width=.7) + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		coord_flip() + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="", fill = "somy") + 
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 8, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(), legend.position = "none",
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank())
		i = i + 2
	}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

  leg1 <- g_legend(tmplot1)
  leg2 <- g_legend(tmplot2)


xp1 = 0.25
yp1 = 0.31
yp2 = 0.325
wc1 = 0.18
diffp = 0.03
wc2 = 0.03
hc1 = 0.30
hc2 = 0.25
hc3 = 0.29

	postscript("~/plotsave/new_fig4_with_everything.eps", width = 8, height = 7, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")

	ggdraw() +
	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +
	draw_plot(plots[[3]],  0, 	yp1*2, wc1, hc1) +
	draw_plot(plots[[7]], xp1, 	yp1*2, wc1, hc1) +
	draw_plot(plots[[9]], xp1*2, yp1*2, wc1, hc1) +
	draw_plot(plots[[11]],xp1*3, yp1*2, wc1, hc1) +

	draw_plot(plots[[15]], 0, 	yp1, wc1, hc1) +
	draw_plot(plots[[17]], xp1,	yp1, wc1, hc1) +
	draw_plot(plots[[19]], xp1*2,	yp1, wc1, hc1) +
	draw_plot(plots[[21]], xp1*3,	yp1, wc1, hc1) +

	draw_plot(plots[[4]], xp1-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +	
	draw_plot(plots[[8]], (xp1*2)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +
	draw_plot(plots[[10]], (xp1*3)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +
	draw_plot(plots[[12]], (xp1*4)-diffp, (yp2*2)-0.005, wc2, hc2-0.005) +

	draw_plot(plots[[16]], xp1-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[18]], (xp1*2)-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[20]], (xp1*3)-diffp, yp2+0.01, wc2, hc2-0.007) +
	draw_plot(plots[[22]], (xp1*4)-diffp, yp2+0.01, wc2, hc2-0.007) +

	draw_plot(f4b, 0, 0, 1, 0.26) +
  	draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.26), size = 16)

    dev.off()
