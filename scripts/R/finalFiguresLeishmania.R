#########################################################
# 						Libraries 						#
#########################################################

library(ggplot2)
library(reshape)
library(plyr)
library(gtools)
library(gdata)
library(grDevices)
library(ggbiplot)
library(Rmisc)
library(grid)


#########################################################
# 						Plotting config 				#
#########################################################


science_theme_jitter = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"), legend.position = c(0.85, 
        0.7), text = element_text(size = 14), axis.text.x = element_text(angle=90))
science_theme3 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"),text = element_text(size = 14), legend.direction = 'horizontal', legend.position = 'top')
science_theme2 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), panel.grid.major = element_line(color = "white"),
    axis.line = element_line(size = 0.7, color = "black"),text = element_text(size = 14))



##################################################################################################################################
#						 	Panel 1 							#
##################################################################################################################################

#######################
# Fig1-A
#######################

# Fuckign correct the compressed chromosome 34 used for normalization #
load("~/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/redefine_piles/.RData")
test.csv = read.csv("~/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/piles/coverage_matrix.txt")

for(i in 1:dim(test.csv)[1]) { #forach sample get 5 llowest
	sample = array(test.csv[i,1])
	med = median(sort(as.numeric(test.csv[i,2:38]))[1:5])
	newval = test.csv[i,"Ld34"]
	coverage[sample,"Ld34"] <- newval/med
}
coverage.melted2 = melt(coverage)
coverage.melted2$variable = factor(coverage.melted2$variable, levels = mixedsort(levels(coverage.melted2$variable)), ordered = TRUE)
newchrs <- gsub(pattern="Ld", coverage.melted2$variable, replacement = "")
coverage.melted2$variable <- factor(newchrs, levels = mixedsort(unique(newchrs)), ordered = TRUE)

# Now fucking plot again
svg(filename = "~/plotsave/coverage_jitter_isolates.svg", width = 8, height = 4)
ggplot(subset(coverage.melted2, variable %in% levels(coverage.melted2$variable)[1:36]), aes( x = variable, y = value*2)) + 
geom_jitter(size = .7) + 
labs(title="", x = "Chromosome", y = "Somy") + theme_bw(base_size = 12, base_family = "Helvetica") + science_theme_jitter
dev.off()

#######################
# Fig1-A
#######################
# Read 

# Transform
coverage.class = gsub(row.names(coverage), pattern = "(^[^0-9]+).*$", replacement = "\\1")
coverage.pca <- prcomp(coverage[1:36], scale. = FALSE)

# Plot
svg(filename = "~/plotsave/pca_biplot_isolates.svg", width = 6, height = 4)
ggbiplot(coverage.pca, labels = NULL, obs.scale = 1, var.scale = 1, var.axes = FALSE,
  groups = coverage.class, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  science_theme3
dev.off()


#################################################################
#						 	Panel 3 							#
#################################################################

# Passages Boxplot

# Passages Density plots - CM: Added chromosome 9
test.p = read.table("~/projects/kinetoplastids/Ldonovani/analysis/SNP_new/DNA/new_passages/Hp0_Hp10_Hp21_intersect_LDBPK.bed")

test.p = read.table("~/projects/kinetoplastids/Ldonovani/analysis/SNP/new_passages_smiley/fixed_new_merged.txt")

colnames(test.p) <- c("Chromosome", "Position", "Frequency", "Sample", "Depth")
chrv <- array(test.p$Chromosome)
test.p$Chromosome <- factor(chrv, levels=mixedsort(unique(chrv)), ordered = TRUE)
test.p$Sample <- array(test.p$Sample)
test.p[which(test.p$Sample == "Ht0"), "Sample"] <- "Splenic amastigote"
test.p[which(test.p$Sample == "Hp0"), "Sample"] <- "P2"
test.p[which(test.p$Sample == "Hp10"), "Sample"] <- "P10"
test.p[which(test.p$Sample == "Hp21"), "Sample"] <- "P20"
test.p$Sample <- factor(test.p$Sample, levels = c("Splenic amastigote", "P2", "P10", "P20"), ordered = TRUE)
svg(filename="~/plotsave/passage_density.svg", height = 4, width = 16)
ggplot(subset(test.p, Chromosome %in% c("5", "9", "14", "15", "20", "26", "31")), aes(x = Frequency, color = Sample)) + geom_density() + facet_wrap(~Chromosome, ncol=7) + labs(x = "Frequency", y = "Density", title="") + theme_bw(base_size = 12, base_family = "Helvetica") + science_theme2
dev.off()

svg(filename="~/plotsave/passage_density_full.svg", height = 20, width = 7)
ggplot(subset(test.p, !Chromosome %in% c("MC", "MI")), aes(x = Frequency)) + geom_density() + facet_grid(Chromosome~Sample) + labs(x = "Frequency", y = "Density", title="") + theme_bw(base_size = 12, base_family = "Helvetica")
dev.off()

# Haplotype map
df.sites = read.table("final.merged.txt")
#df.sites = read.table("merged.for20.txt2")
# needed to do`§
df.sites = read.table("merged.s.txt")
df.sites = read.table("~/Desktop/merged.wp20.grey.txt")
colnames(df.sites) = c("chr", "position", "CL10", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "P20")
df.sites[which(df.sites$chr == "3"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9", "P20")] <- 0
df.sites[which(df.sites$chr == "8"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9")] <- 0
df.sites[which(df.sites$chr == "14"),c("CL1", "CL8")] <- 0
df.sites[which(df.sites$chr == "15"),c("CL8")] <- 0
df.sites[which(df.sites$chr == "16"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9")] <- 0
df.sites[which(df.sites$chr == "20"),c("CL1", "CL8")] <- 0
df.sites[which(df.sites$chr == "23"),c("CL1", "CL8")] <- 0
df.sites[which(df.sites$chr == "32"),c("CL10", "CL3", "CL4", "CL6", "CL7", "CL9", "P20")] <- 0


pdf("~/plotsave/haplotype_map_clones.pdf")
for(c in chromosomes) {
	df.sub = subset(df.sites, chr == c)
	df.sub$npos  = seq(from=1, to=dim(df.sub)[1])
	df.sites.melt = melt(df.sub, id.vars=c("chr", "position", "npos"))
	print(ggplot(df.sites.melt, aes(variable, npos)) + geom_tile(aes(fill=factor(value))) + scale_fill_manual(values=colors) + labs(title = c, x = "Clone", y = "Relative position site", fill = "Allele") + coord_flip() )
}
dev.off()


svg(filename="~/plotsave/haplotype_map_clones_global2.svg", height = 4, width = 14)
#print(ggplot(subset(df.sites.melt, chr %in% c("20")), aes(y=factor(variable), x=npos)) + 
print(ggplot(subset(df.sites.melt, chr %in% c("5", "9", "14", "15","26", "31")), aes(y=factor(variable), x=npos)) + 
	geom_tile(aes(fill=factor(value), height = .7)) + 
	scale_fill_manual(values=colors, labels=c("A","C","G","T"), breaks = c(0,1,2,3)) +
	scale_x_continuous(expand = c(0, 0)) + 
	facet_grid(.~chr, scales = "free") +
	labs(title = "Genome haplotype map", y = "Sample", x = "Genome position", fill = "Major allele") +
	theme_bw(base_size = 12, base_family = "Helvetica") +
	theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
	#theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
	#science_theme2
	)
dev.off()
pdf("~/plotsave/haplotype_map_clones_global.pdf", height = 4, width = 10)
for (c in levels(df.sites.melt$chr)) {
print(ggplot(subset(df.sites.melt, chr == c), aes(y=factor(variable), x=npos)) + 
	geom_tile(aes(fill=factor(value), height = .7)) + 
	scale_fill_manual(values=colors, labels=c("A","C","G","T"), breaks = c(0:3)) +
	scale_x_continuous(expand = c(0, 0)) + 
	facet_grid(.~chr, scales = "free") +
	labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
	theme_bw(base_size = 12, base_family = "Helvetica") +
	theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
	#theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
	#science_theme2
	)
}
dev.off()


### Prepare merged piles from p21 and clones
df.passages.piles <- read.table("~/Hp0_Hp10_Hp21_intersect_LDBPK.bed")
colnames(df.passages.piles) <- c("chromosome", "start", "end", "sample","frequency")

df.clones.piles <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/analysis/SNP_new/DNA/DNAseq_2014/pileup_smiley/DNA_clones.Rpile.txt")
df.clones.piles <- read.table("DNA_clones.Rpile.txt")



# Almost the last figure: Hapalotype with
svg("~/plotsave/clones_densities_test.svg", height = 10, width = 2)
ggplot(subset(df.clones.piles, chromosome == "31"), aes(frequency)) + geom_density() + facet_wrap(~sample, ncol = 1 ) + labs(x = "", y="")
dev.off()


df.t <- subset(df.passages.piles, !(chromosome %in% c("MC", "MI")) & sample == "Hp21")[,c("chromosome", "sample", "frequency")]
df.t$sample <- gsub(df.t$sample, pattern = "H", replacement = "")
df.t2 <- subset(df.clones.piles, !(chromosome %in% c("MC", "MI")) & !(sample %in% c("DNA_0", "DNA_AX")))[,c("chromosome", "sample", "frequency")]
df.t2$sample <- gsub(df.t2$sample, pattern = "DNA_", replacement = "")
df.pclones.piles <- rbind(df.t, df.t2)

# Now merge with somy numbers
colnames(gpkm.clones.m.medians) <- c("sample", "chromosome", "somy")
df.pclones.piles.medians <- merge(subset(gpkm.clones.m.medians, !(sample %in% c("Ax_ama", "Sp_ama", "OSp_ama")) & !(chromosome %in% c("MC", "MI"))), df.pclones.piles)
newvar = array(df.pclones.piles.medians$sample)
df.pclones.piles.medians$sample <- factor(newvar, levels = rev(c("p21","CL1", "CL8", "CL4", "CL3", "CL6", "CL7","CL9", "CL10")), ordered = TRUE)

# 
df.clones.piles.density <- subset(df.clones.piles, !(chromosome %in% c("MC", "MI")) & !(sample %in% c("DNA_0", "DNA_AX")))[,c("chromosome", "sample", "frequency")]
df.clones.piles.density$sample <- gsub(df.t2$sample, pattern = "DNA_", replacement = "")
newvar = array(df.clones.piles.density$sample)
df.clones.piles.density$sample <- factor(newvar, levels = rev(c("CL1", "CL8", "CL4", "CL3", "CL6", "CL7","CL9", "CL10")), ordered = TRUE)

chrv <- array(df.clones.piles.density$chromosome)
df.clones.piles.density$chromosome <- factor(chrv, levels=mixedsort(unique(chrv)), ordered = TRUE)

svg(filename="~/plotsave/clones_density.svg", height = 20, width = 10)
ggplot(df.clones.piles.density, aes(x = frequency)) + geom_density() + facet_grid(chromosome~sample) + labs(x = "Frequency", y = "Density", title="") + theme_bw(base_size = 12, base_family = "Helvetica") 
dev.off()


# Testing join haplo+density plot
pdf("~/plotsave/haplotype_density_combo.pdf", height = 5, width = 6)
for (c in levels(df.sites.melt$chr)) {

	plots <- list()  
	layout <- matrix(c(1, 1, 2), nrow = 1, byrow = TRUE)

	plots[[1]] <- ggplot(subset(df.sites.melt, chr == c), aes(y=factor(variable), x=npos)) + 
	geom_tile(aes(fill=factor(value), height = .7)) + 
	scale_fill_manual(values=colors, labels=c("A","C","G","T"), breaks = c(0:3)) +
	scale_x_continuous(expand = c(0, 0)) + 
	facet_grid(.~chr, scales = "free") +
	#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
	labs(title = "", y = "", x = "", fill = "") +
	theme_bw(base_size = 12, base_family = "Helvetica") +
	theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")), legend.position = "none")
	#theme(panel.grid=element_blank(), panel.border=element_blank(), axis.text = element_text(margin = unit(2,"lines")))
	#science_theme2
	plots[[2]] <- ggplot(subset(df.pclones.piles.medians, chromosome == c), aes(frequency, fill = somy)) + 
	geom_density() + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
	facet_wrap(~sample, ncol = 1 ) + 
	#labs(x = "", y="", fill = "Somy") + 
	labs(x = "", y="", fill = "") + 
	theme(plot.margin = unit(c(1.8,1,0.2,0), "cm"), panel.margin.y = unit(0.4, "lines"), 
		axis.text.x = element_text(angle = 90), strip.background = element_blank(), strip.text.x = element_blank(), legend.position = "none")
	multiplot(plotlist = plots, layout = layout)
}
dev.off()

# Missing data from ploidy gpkm median values for every sample
# TODO: After get another table of relative fuckind depth and plot side by side too

#################################################################
#						 	Panel 3 							#
#################################################################

# Add liver vs spleen boxplot in the supplement

# Here it has be colored by aneuploidy frequency?
# Dotplot RNA vs DNA 
melted.rpkm.t = melt(rpkm[,c(2,8:15)], id.vars=c("CHROMOSOME"))
melted.gpkm.t = melt(gpkm[,c(2,8:15)], id.vars=c("CHROMOSOME"))
medians.melted.t = ddply(melted.rpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})[,1:2]
medians.melted.t$rpkm = ddply(melted.rpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})$V1
medians.melted.t$gpkm = ddply(melted.gpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})$V1
chrv <- array(medians.melted.t$CHROMOSOME)
medians.melted.t$CHROMOSOME <- factor(chrv, levels=mixedsort(unique(chrv)))
svg(filename="~/plotsave/median_dna_rna_clones.svg", height = 4, width = 6)
ggplot(subset(medians.melted.t, !(CHROMOSOME %in% c("MC", "MI"))), aes(gpkm, rpkm)) + 
geom_point() + 
geom_point() +
geom_abline(slope=1, intercept=-30) + geom_point(data=subset(medians.melted.t, (CHROMOSOME %in% c("5", "14", "15", "20", "26", "31"))), aes(color=CHROMOSOME)) + labs(color = "Chromosome", x="Genomic coverage", y="Transcriptomic coverage", title="") + theme_bw(base_size = 12, base_family = "Helvetica") + science_theme2
dev.off()

# Plot for SNP numbers on , coloring the dots using frequency of aneupoidy


read.table("samples_to_snp_count.counts.txt") -> df.snps.field
colnames(df.snps.field) <- c("chromosome", "sites", "length")
vars <- gsub(df.snps.field$chromosome, pattern="\\|\\S+", replacement="")
df.snps.field$chromosome <- factor(vars, levels = unique(mixedsort(vars)), ordered = TRUE)
num_samples <- length(levels(coverage.melted$Isolate))
df.aneuploid.field <- ddply(coverage.melted, c("variable"), function(x){ length(which(x$value > 1.2)) / num_samples })[1:36,]
chrreplace <- gsub(df.aneuploid.field$variable, pattern="Ld", replacement="")
df.aneuploid.field$variable <- factor(chrreplace, levels = , ordered = TRUE)

df.snps.field$anefreq <- rep(0,38)
df.snps.field[3:38,"anefreq"] <- df.aneuploid.field$V1
science_theme2 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"),text = element_text(size = 12))
svg('~/plotsave/snp_count_isolates2.svg', height = 6, width = 10)
ggplot(df.snps.field, aes(length/1000000, sites, label=chromosome)) + geom_point(size = 6, shape = 16, aes(colour = anefreq)) + geom_text(size=3) + 
geom_abline(color="blue", slope = 0.0004*1000000) + geom_abline(color="green", slope = 0.0007*1000000) + geom_abline(color="red", slope = 0.0018*1000000) +
scale_x_continuous(breaks=c(1,2),labels=c("1","2")) +
#scale_shape_discrete(solid=T) +
scale_colour_gradient(low = "white", high = "#EE3B3B") +
labs(colour = "Aneuploidy frequency", x="Chromosome length (Mb)", y="# Het variant sites") + 
theme_bw(base_size = 12, base_family = "Helvetica") + science_theme2
dev.off()

#################################################################



# Other crappy staff
# Different haplotypes in liver and spleen
read.table("C7A9MANXX_1_15nf.filtered.bed") -> df.test.1
read.table("C7A9MANXX_4_8nf.filtered.bed") -> df.test.4
pdf('test4.pdf')
ggplot(subset(df.test.4, V1 == "20|FR799607.2" & V3 > 0.4 & V3 < 0.5),aes(V3)) + geom_histogram() + coord_cartesian(xlim=c(0,1), ylim=c(0,5))
dev.off()
pos = unique(subset(df.test.4, V1 == "20|FR799607.2" & V3 > 0.4 & V3 < 0.5)$V2)
pdf('test1.pdf')
ggplot(subset(df.test.1, V1 == "20|FR799607.2" & V2 %in% pos),aes(V3)) + geom_histogram() + coord_cartesian(xlim=c(0,1), ylim=c(0,5))
dev.off()


# Full profiles
df.test.1$sample = "spleen"
df.test.4$sample = "liver"
df.spleen.liver <- rbind(df.test.1, df.test.4)
df.spleen.liver$V1<- gsub(df.spleen.liver$V1, pattern="\\|\\S+", replace="")
df.spleen.liver$V1 <- factor( df.spleen.liver$V1, levels = mixedsort(unique(array(df.spleen.liver$V1))), ordered = TRUE )
colnames(df.spleen.liver) <- c("chromosome", "position", "frequency", "sample")
svg(filename="~/plotsave/liver_spleen_all_density.svg", height = 20, width = 4)
ggplot(subset(df.spleen.liver, !chromosome %in% c("MC", "MI")),aes(frequency)) + geom_density() + facet_grid(chromosome~sample, scales="free")  + coord_cartesian(xlim=c(0,1)) + labs(x = "Frequency", y = "Density", title="") + theme_bw(base_size = 12, base_family = "Helvetica") + science_theme2
dev.off()
####
reads2rpkm <- function (m_reads, tr_length, total_reads) {
  
  RPK <- (m_reads / (tr_length / 1000)) # mapped reads / transcript length in kb
  RPKM <- (RPK / (total_reads /  1000000)) # RPK / total num reads in million
  
  return (RPKM)
}
df.gene_counts.DNA.passages <- read.table("/Users/admin/work/Ldonovani/passage_number/axam_Ht0_Hp0_Hp10_Hp21.bed", header=T)
df.gene_counts.DNA.passages$Chromsome = gsub(df.gene_counts.DNA.passages$Chromsome, pattern="(.*)\\|F.*", replacement="\\1")
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
svg(filename="~/plotsave/passage_boxplot.svg", height = 4, width = 14)
ggplot(subset(gpkm.passages.m, Chromosome != "MC" & variable %in% c("Sp_ama", "p2", "p10", "p20")), aes(Chromosome, value, fill = variable)) + 
geom_boxplot( outlier.shape = NA) + facet_grid(.~Chromosome, scales = "free_x") + coord_cartesian(ylim=c(0,160))  + 
scale_fill_manual(values = gg_color_hue(4),labels = c("Splenic amastigote", "P2", "P10", "P20")) +
labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + science_theme2
dev.off()

# Now let's do the same with Clones
gpkm.clones <- gpkm
colnames(gpkm.clones) <- c("Gene","Chromosome","Start","End","Strand","Ax_ama","Sp_ama","CL1" ,"CL3","CL4","CL6","CL7", "CL8", "CL9", "CL10", "OSp_ama")
gpkm.clones$Chromosome = gsub(gpkm.clones$Chromosome, pattern="(.*)\\|F.*", replacement="\\1")
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

svg(filename="~/plotsave/clones_subset_boxplot.svg", height = 4, width = 10)
ggplot(subset(gpkm.clones.m, Chromosome %in% c("5", "9", "14", "15", "20", "26", "31") & variable %in% c("p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")), aes(Chromosome, value, fill = variable)) + geom_boxplot( outlier.shape = NA) + coord_cartesian(ylim=c(0,160)) + facet_grid(.~Chromosome, scales = "free_x") + labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
scale_fill_manual(values = gg_color_hue(9),labels = c("P20", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")) +
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + science_theme2
dev.off()

svg(filename="~/plotsave/clones_boxplot.svg", height = 4, width = 14)
ggplot(subset(gpkm.clones.m, variable %in% c("CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")), aes(Chromosome, value, fill = variable)) + geom_boxplot( outlier.shape = NA, coef = 0 ) + 
coord_cartesian(ylim=c(0,160)) + facet_grid(.~Chromosome, scales = "free_x") + 
labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + science_theme2
dev.off()

rpkm.clones <- rpkm
colnames(rpkm.clones) <- c("Gene","Chromosome","Start","End","Strand","Ax_ama","Sp_ama","CL1" ,"CL3","CL4","CL6","CL7", "CL8", "CL9", "CL10", "OSp_ama")
rpkm.clones$Chromosome = gsub(rpkm.clones$Chromosome, pattern="(.*)\\|F.*", replacement="\\1")
rpkm.clones$Chromosome <- factor(rpkm.clones$Chromosome, levels = mixedsort(unique(rpkm.clones$Chromosome)), ordered = TRUE)
rpkm.clones.m <- melt(rpkm.clones, measure.vars = c("Ax_ama", "Sp_ama", "OSp_ama", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10"))

svg(filename="~/plotsave/clones_rna_boxplot.svg", height = 4, width = 14)
ggplot(subset(rpkm.clones.m, variable %in% c("CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10") & !(Chromosome %in% c("MC", "MI"))), aes(Chromosome, value, fill = variable)) + geom_boxplot( outlier.shape = NA, coef = 0 ) + 
coord_cartesian(ylim=c(0,160)) + facet_grid(.~Chromosome, scales = "free_x") + 
labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + science_theme2
dev.off()


#################################################################################################################################
#														Supplementary															#
#################################################################################################################################

df.gene_counts.DNA.spleen_liver <- read.table("~/Desktop/C7A9MANXX_1_15nf.C7A9MANXX_4_8nf.bed", header=T)
df.gene_counts.DNA.spleen_liver$Chromosome = gsub(df.gene_counts.DNA.spleen_liver$Chromosome, pattern="(.*)\\|F.*", replacement="\\1")
colnames(df.gene_counts.DNA.spleen_liver) <- c("Geneid","Chromosome","Start","End","Score","Strand","Spleen","Liver")
gpkm.spleen_liver <- df.gene_counts.DNA.spleen_liver
gpkm.spleen_liver$Chromosome <- factor(gpkm.spleen_liver$Chromosome, levels = mixedsort(unique(gpkm.spleen_liver$Chromosome)), ordered = TRUE)

T <- sum(gpkm.spleen_liver$Spleen)
gpkm.spleen_liver$Spleen<-apply(gpkm.spleen_liver,1, function(x){ reads2rpkm(as.numeric(x[7]),(as.numeric(x[4])-as.numeric(x[3])),sum(gpkm.spleen_liver$Spleen))})
T <- sum(gpkm.spleen_liver$Liver)
gpkm.spleen_liver$Liver<-apply(gpkm.spleen_liver,1, function(x){ reads2rpkm(as.numeric(x[8]),(as.numeric(x[4])-as.numeric(x[3])),sum(gpkm.spleen_liver$Liver))})

gpkm.spleen_liver.m <- melt(gpkm.spleen_liver, measure.vars = c("Spleen", "Liver"))
svg(filename="~/plotsave/spleen_liver_boxplot.svg", height = 4, width = 14)
ggplot(subset(gpkm.spleen_liver.m, Chromosome != "MC"), aes(Chromosome, value, fill = variable)) + geom_boxplot( outlier.shape = NA) + facet_grid(.~Chromosome, scales = "free_x") + coord_cartesian(ylim=c(0,160))  + labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) + science_theme2
dev.off()



############# Suplementary figure as for the deep shit histograms of diploid/trisomic/etc ...

#################################################################################################################################
#														Stacked disomic/trisomic profiles field isolates															#
#################################################################################################################################

	df.het <- read.table("all.pile.Rpile.depth.1.6.smiley.new.txt",header=F)
	colnames(df.het) <- c("chromosome", "start","frequency","name", "depth")
	newchrs <- gsub(pattern="Ld", df.het$chromosome, replacement = "")
	df.het$chromosome = factor( newchrs, levels= mixedsort(unique(newchrs)), ordered=TRUE )

	###
	# What if you want to facet all of them?
	# With all of them?

	df.het.by.somy <- ddply(df.het, .(chromosome), function(d) {
		chr = unique(array(d$chromosome))
		if(chr %in% c("MC", "MI")) {
			return()
		}
		array_coverage 		= get(paste("Ld", chr, sep = ""), coverage)
		array_isolates 		= array(rownames(coverage))

		monosomic_samples   = array_isolates[array_coverage < 0.85 & !is.na(array_coverage)]
		disomic_samples     = array_isolates[array_coverage > 0.85 & array_coverage < 1.15 & !is.na(array_coverage)]
		trisomic_samples 	= array_isolates[array_coverage > 1.15 & array_coverage < 1.85 & !is.na(array_coverage)]
		tetrasomic_samples 	= array_isolates[array_coverage > 1.85 & array_coverage < 2.25 & !is.na(array_coverage)]
		pentasomic_samples  = array_isolates[array_coverage > 2.25 & !is.na(array_coverage)]

		somy_samples <- {}

		if(length(monosomic_samples)) {
			somy_samples <- rbind(somy_samples,  cbind(subset(d, name %in% monosomic_samples), somy = "monosomic"))
		}
		if(length(disomic_samples)) {
			somy_samples <- rbind(somy_samples,  cbind(subset(d, name %in% disomic_samples), somy = "disomic"))
		}
		if(length(trisomic_samples)) {
			somy_samples <- rbind(somy_samples,  cbind(subset(d, name %in% trisomic_samples), somy = "trisomic"))
		}
		if(length(tetrasomic_samples)) {
			somy_samples <- rbind(somy_samples,  cbind(subset(d, name %in% tetrasomic_samples), somy = "tetrasomic"))
		}
		if(length(pentasomic_samples)) {
			somy_samples <- rbind(somy_samples,  cbind(subset(d, name %in% pentasomic_samples), somy = "pentasomic"))
		}
		if( length(monosomic_samples) | length(disomic_samples) | length(trisomic_samples) | length(tetrasomic_samples) | length(pentasomic_samples)) {
			return(somy_samples)
		} else {
			return
		}
	})

	count_samples_bins <- function (df.ss) {
				df.new <- df.ss
				df.new$cuts <- cut(df.new$frequency, breaks = seq(0, 1, by =0.1))
				samples <- ddply(df.new, .(cuts, somy,chromosome), function(df){length(unique(df$name))})
				colnames(samples)[3] <- "samples"
				return(samples)
			}

	count_samples <- function (df.ss) {
		df.new <- df.ss
		samples <- ddply(df.new, .(somy, chromosome), function(df){length(unique(df$name))})
		colnames(samples)[3] <- "samples"
		return(samples)
	}

	df.somy.counts <- count_samples(df.het.by.somy)

	df.somy.counts$somy <- factor(df.somy.counts$somy, levels = c("none", "monosomic", "disomic", "trisomic", "tetrasomic", "pentasomic"), ordered = TRUE)
	for(c in levels(df.somy.counts$chromosome)) {
		print(c)
		tot = 204 - sum(as.integer(subset(df.somy.counts, chromosome == c)$samples))
		df.somy.counts <- rbind(df.somy.counts, c("nonusable", c, tot))
	}
	df.somy.counts$samples <- as.integer(df.somy.counts$samples)
	df.somy.counts <- subset(df.somy.counts, !chromosome %in% c("MC", "MI"))
	df.somy.counts$chromosome <- factor(df.somy.counts$chromosome, unique(mixedsort(array(df.somy.counts$chromosome))), ordered = TRUE)

	svg(filename = "~/plotsave/chr_somy_samples.svg", width = 12, height = 8)	
	ggplot(df.somy.counts, aes(chromosome, samples, fill = somy)) + geom_bar(stat = "identity") +
	labs(x= "Chromosome", y = "Number of samples", fill = "Somy") +
	theme_bw(base_size = 12, base_family = "Helvetica") +
	science_theme2
	dev.off()

	svg(filename = "~/plotsave/chr_disomic_histogram.svg", width = 8, height = 20)
		ggplot(subset(df.het.by.somy, somy %in% c("any", "disomic", "trisomic")), aes(frequency)) + 
		geom_histogram(alpha=.5,binwidth=0.01) + 
		facet_grid(chromosome~somy, scales="free_y") + 
		theme(axis.text.x = element_text(angle=90)) +
		theme_bw()
	dev.off()

	svg(filename = "~/plotsave/chr_disomic_histogram_sample_counts.svg", width = 8, height = 4)
		ggplot(df.disomic.counts, aes(cuts, samples)) + geom_bar(stat="identity") + 
		coord_cartesian(ylim = c(0,200)) + 
		facet_wrap(~chromosome, scales="free", ncol=3) + 
		labs(x="frequency", y = "# Samples") + 
		theme(axis.text.x = element_text(angle=90))
	dev.off()

	svg(filename = "~/plotsave/chr_trisomic_histogram.svg", width = 8, height = 4)
		ggplot(subet(df.het.by.somy, somy == "trisomic"), aes(frequency)) + 
		geom_histogram(alpha=.5,binwidth=0.01) + 
		facet_wrap(~chromosome, scales="free", ncol=3) + 
		theme(axis.text.x = element_text(angle=90))
	dev.off()

	svg(filename = "~/plotsave/chr_trisomic_histogram_sample_counts.svg", width = 8, height = 4)
		ggplot(df.trisomic.counts, aes(cuts, samples)) + geom_bar(stat="identity") + 
		coord_cartesian(ylim = c(0,200)) + 
		facet_wrap(~chromosome, scales="free", ncol=3) + 
		labs(x="frequency", y = "# Samples") + 
		theme(axis.text.x = element_text(angle=90))
	dev.off()


### Get a couple of samples with profiles of trisomic samples in chromosomes 5,11,12
union_5_11 = intersect( 	array(unique(subset(df.het.by.somy, chromosome == 5 & somy == "trisomic")$name)), 
		array(unique(subset(df.het.by.somy, chromosome == 11 & somy == "trisomic")$name)))

union_5_11_12 = intersect(	union_5_11, array(unique(subset(df.het.by.somy, chromosome == 12 & somy == "trisomic")$name)) )

svg("~/tmp/chr_5_11_12_isolate_frequency_histogram.svg", width = 14, height = 6)
ggplot(subset(df.het.by.somy, chromosome %in% c("5","11","12") & name %in% union_5_11_12), aes(frequency)) +
+ geom_histogram() + geom_density() + facet_grid(chromosome~name, scales = "free_y") +
+ theme_bw(base_size = 12, base_family = "Helvetica") + science_theme2
dev.off()