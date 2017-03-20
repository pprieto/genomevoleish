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
library(corrplot)
library(gridExtra)

# Theme stuff
science_theme_jitter = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"), legend.position = c(0.85, 
        0.7), axis.text.x = element_text(angle=90))
science_theme3 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"), legend.direction = 'horizontal', legend.position = 'top')
science_theme2 = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), panel.grid.major = element_line(color = "white"),
    axis.line = element_line(size = 0.7, color = "black"))

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

f1a <- ggplot(subset(coverage.melted2, variable %in% levels(coverage.melted2$variable)[1:36]), aes( x = variable, y = value*2)) + 
	geom_jitter(size = .7) + 
	labs(title="", x = "Chromosome", y = "Somy") + theme_bw(base_size = 8, base_family = "ArialMT") #+ science_theme_jitter

number_ticks <- function(n) {function(limits) pretty(limits, n)}
	
new_f1a <- ggplot(subset(coverage.melted2, variable %in% levels(coverage.melted2$variable)[1:36]), aes( x = value*2)) + 
	geom_density() +
	facet_grid(variable~., scales = "free") +
	scale_y_continuous(breaks=number_ticks(2)) +	
	labs(title="", x = "Somy", y = "Frequency") + 
	theme_classic(base_size = 12, base_family = "ArialMT") +
	theme(text = element_text(face = "bold"), strip.text.y = element_text(angle = 0), strip.background = element_blank(), axis.text = element_text(size = 8))

coverage = read.csv("~pprieto/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/piles/coverage_matrix.txt", header = TRUE)
# Normalize the coverages
for(i in 1:dim(coverage)[1]) {
	im = median(array(unlist(coverage[i,2:38])), na.rm = T)
	for(j in 2:38) {
		coverage[i,j] = (coverage[i,j]/im)
	}
}

row.names(coverage) <- coverage$Isolate
coverage$Isolate <- {}

coverage.class = gsub(row.names(coverage), pattern = "(^[^0-9]+).*$", replacement = "\\1")
coverage.pca <- prcomp(coverage[1:36], scale. = FALSE)

f1b <- ggbiplot::ggbiplot(coverage.pca, labels = NULL, obs.scale = 0, var.scale = 1, var.axes = FALSE,
  groups = coverage.class, ellipse = TRUE, circle = TRUE) +#+
#  scale_color_discrete(name = '') +
	theme_bw(base_size = 8, base_family = "ArialMT") +
  science_theme3


postscript("~/tmp/fig1.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
ggdraw() +
	draw_plot(f1a, 0, .5, 1, .5) +
  	draw_plot(f1b, 0, 0, .6, .5) +
  	draw_plot(f1c1, 0.5, 0, (0.50/3)*2, .5) +
  	draw_plot(f1c2, 0.5 + (0.50/3)*2, 0, (0.50/3), .48) +
  	draw_plot_label(c("a", "b", "c"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)
dev.off()


#fig1c <- function() {
	df.het <- read.table("~/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/redefine_piles/all.pile.Rpile.depth.1.6.smiley.new.txt",header=F)
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

		monosomic_samples   = array_isolates[array_coverage < 0.85]
		disomic_samples     = array_isolates[array_coverage > 0.85 & array_coverage < 1.15]
		trisomic_samples 	= array_isolates[array_coverage > 1.15 & array_coverage < 1.85]
		tetrasomic_samples 	= array_isolates[array_coverage > 1.85 & array_coverage < 2.25]
		pentasomic_samples  = array_isolates[array_coverage > 2.25]

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

	df.allele.std <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/contrast_alleles/all/all_real_results.txt")
	colnames(df.allele.std) <- c("chromosome", "start", "allele", "frequency")
	df.allele.std$somy <- "std"
	df.het.by.somy.plot <- rbind(df.het.by.somy[,c("chromosome", "start", "frequency", "somy")], df.allele.std[,c("chromosome", "start", "frequency", "somy")])

	
	df.het.by.somy.plot.vline <- expand.grid(levels(df.het.by.somy.plot$chromosome) ,levels(df.het.by.somy.plot$somy))
	colnames(df.het.by.somy.plot.vline) <- c("chromosome", "somy")
	df.het.by.somy.plot.vline$hline <- c(rep(0.5, 76*2), rep(0.07, 76))

	
	f1c1 <- ggplot(subset(df.het.by.somy.plot, chromosome %in% c("5", "11", "12") & somy %in% c("disomic", "trisomic")), aes(frequency)) + 
		geom_density(adjust=1.5) +
		geom_vline(data = subset(df.het.by.somy.plot.vline, chromosome %in% c("5", "11", "12") & somy %in% c("disomic", "trisomic")), aes(xintercept = hline), color = "red") + 
		facet_grid(chromosome~somy, scales="free") + 
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(strip.text.y = element_blank(), strip.background = element_blank())

	f1c2 <- ggplot(subset(df.het.by.somy.plot, chromosome %in% c("5", "11", "12") & somy %in% c("std")), aes(frequency)) + 
		geom_histogram(binwidth=0.01) + 
		geom_vline(data = subset(df.het.by.somy.plot.vline, chromosome %in% c("5", "11", "12") & somy %in% c("std")), aes(xintercept = hline), color = "red") + 
		scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
		facet_grid(chromosome~somy, scales="free_y") + 
		labs(x = "std") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())

		#### New STD plot for all the chromosomes
		if(FALSE) {
			f1c1 <- ggplot(subset(df.het.by.somy.plot, !(chromosome %in% c("MC", "MI")) & somy %in% c("disomic", "trisomic")), aes(frequency)) + 
				geom_density(adjust=1.5) +
				geom_vline(data = subset(df.het.by.somy.plot.vline, !(chromosome %in% c("MC", "MI")) & somy %in% c("disomic", "trisomic")), aes(xintercept = hline), color = "red") + 
				facet_grid(chromosome~somy, scales="free") + 
				theme_bw(base_size = 8, base_family = "ArialMT") +
				theme(strip.text.y = element_blank(), strip.background = element_blank())

			f1c2 <- ggplot(subset(df.het.by.somy.plot, !(chromosome %in% c("MC", "MI")) & somy %in% c("std")), aes(frequency)) + 
				geom_histogram(binwidth=0.01) + 
				geom_vline(data = subset(df.het.by.somy.plot.vline, !(chromosome %in% c("MC", "MI")) & somy %in% c("std")), aes(xintercept = hline), color = "red") + 
				scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
				facet_grid(chromosome~somy, scales="free_y") + 
				labs(x = "std") +
				theme_bw(base_size = 8, base_family = "ArialMT") +
				theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())

				postscript("~/tmp/newfig1.eps", width = 8, height = 18, horizontal = FALSE, onefile = FALSE, paper = "special", 
				colormodel = "rgb", family = "ArialMT")
			ggdraw() +
			  	draw_plot(f1c1, 0, 0, (1/3)*2, 1) +
			  	draw_plot(f1c2, 0 + (1/3)*2, 0, (1/3), 0.99) +
			  	draw_plot_label(c("c"), c(0), c(0), size = 15)
			dev.off()

			####
			postscript("~/tmp/newfig6.eps", width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special", 
			colormodel = "rgb", family = "ArialMT")

				ggplot(df.snps.field, aes(length/1000000, sites, label=chromosome)) + geom_point(size = 6, shape = 16, aes(colour = anefreq)) + geom_text(size=3) + 
				geom_smooth(data=subset(df.snps.field, anefreq > 0.75 ), method = "lm", colour = "red") +
				geom_smooth(data=subset(df.snps.field, anefreq > 0.25 & anefreq <0.75), method = "lm", colour = "green") +
				geom_smooth(data=subset(df.snps.field, anefreq <0.25), method = "lm", colour = "blue") +
				scale_x_continuous(breaks=c(1,2),labels=c("1","2")) +
				#scale_shape_discrete(solid=T) +
				scale_colour_gradient(low = "white", high = "#EE3B3B") +
				labs(colour = "Aneuploidy\nfrequency", x="Chromosome length (Mb)", y="# Het variant sites") + 
				theme_bw(base_size = 8, base_family = "ArialMT") + science_theme2
			
			dev.off()
 
		}
#	return(f1c)
#}

#func3a <- function() {
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

f3a <- ggplot(subset(gpkm.passages.m, Chromosome %in% c("5", "9", "14", "15", "20", "23", "26", "31") & variable %in% c("Sp_ama", "p2", "p10", "p20")), aes(Chromosome, value, fill = variable)) + 
	geom_boxplot( outlier.shape = NA) + facet_grid(.~Chromosome, scales = "free_x") + coord_cartesian(ylim=c(0,160))  + 
	scale_fill_manual(values = gg_color_hue(4),labels = c("Splenic amastigote", "p2", "p10", "p20")) +
	labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
	science_theme2 +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),  
	     legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') + 
	guides(fill=guide_legend(nrow=1,byrow=FALSE))

f3aS <- ggplot(subset(gpkm.passages.m, Chromosome != "MC" & variable %in% c("Sp_ama", "p2", "p10", "p20")), aes(Chromosome, value, fill = variable)) + 
geom_boxplot( outlier.shape = NA) + facet_grid(.~Chromosome, scales = "free_x") + coord_cartesian(ylim=c(0,160))  + 
scale_fill_manual(values = gg_color_hue(4),labels = c("Splenic amastigote", "P2", "P10", "P20")) +
labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
science_theme2 +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),  
	     legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') + 
	guides(fill=guide_legend(nrow=1,byrow=FALSE))
#return(f3a)
#}


#func3b <- function() {
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

	f3b <- ggplot(subset(gpkm.clones.m, Chromosome %in% c("5", "9", "14", "15", "20", "23", "26", "31") & variable %in% c("p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")), aes(Chromosome, value, fill = variable)) + geom_boxplot( outlier.shape = NA) + coord_cartesian(ylim=c(0,160)) + facet_grid(.~Chromosome, scales = "free_x") + labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
	scale_fill_manual(values = gg_color_hue(9),labels = c("p20", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")) +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	science_theme2 + 
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),  
	     legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') + 
	guides(fill=guide_legend(nrow=1,byrow=FALSE))

#	return(f3b)
#}

postscript("~/plotsave/newfig3.eps", width = 6, height = 5, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
ggdraw() +
	draw_plot(f3a, 0, .5, 1, .5) +
  	 draw_plot(f3b, 0, 0, 1, .5) +
  	draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.5), size = 15)
dev.off()

svg(filename = "~/plotsave/fig3afull.svg", width = 18, height = 12)
	print(f3aS)
dev.off()

#f4a = func4a() {

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

	reads2rpkm <- function (m_reads, tr_length, total_reads) {
	  
	  RPK <- (m_reads / (tr_length / 1000)) # mapped reads / transcript length in kb
	  RPKM <- (RPK / (total_reads /  1000000)) # RPK / total num reads in million
	  
	  return (RPKM)
	}

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

	df.pclones.tiles.count.somy$sample <- factor(df.pclones.tiles.count.somy$sample, levels = rev(levels(df.pclones.tiles.count.somy$sample)), ordered = TRUE)

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
		scale_fill_manual(values=pcolors, labels=c("A","C","G","T", "-")) +
		scale_x_continuous(expand = c(0, 0)) + 
		facet_grid(.~chr, scales = "free") +
		#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
		labs(title = "", y = "", x = "", fill = "alelle") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(strip.background = element_blank(), plot.margin = unit(c(-0.3,0,-0.3,0), "cm"), 
			legend.direction = 'horizontal', legend.position = 'top') + 
			guides(fill=guide_legend(title = "allele", nrow=1,byrow=FALSE))

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


	
	tmplot3 <- ggplot(subset(df.pclones.tiles.count.somy, chromosome == "26"), aes(sample, count, fill = somy)) + 
		geom_bar(stat = "identity", width = .6) + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		coord_flip() + 
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
#}

#f4b = func4b() {
	test.p = read.table("~/projects/kinetoplastids/Ldonovani/analysis/SNP/new_passages_smiley/fixed_new_merged.txt")

	colnames(test.p) <- c("Chromosome", "Position", "Frequency", "Sample", "Depth")
	chrv <- array(test.p$Chromosome)
	test.p$Chromosome <- factor(chrv, levels=mixedsort(unique(chrv)), ordered = TRUE)
	test.p$Sample <- array(test.p$Sample)
	test.p[which(test.p$Sample == "Ht0"), "Sample"] <- "Sp-ama"
	test.p[which(test.p$Sample == "Hp0"), "Sample"] <- "p2"
	test.p[which(test.p$Sample == "Hp10"), "Sample"] <- "p10"
	test.p[which(test.p$Sample == "Hp21"), "Sample"] <- "p20"
	test.p$Sample <- factor(test.p$Sample, levels = c("Sp-ama", "p2", "p10", "p20"), ordered = TRUE)
	
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
#}

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

f5b = func5a() {
	df1 <- data.frame(
		strain = factor(c("p2", "p2", "p10", "p10", "p20", "p20", "LdB", "LdB"), levels = c("p2", "p10", "p20", "LdB"), ordered=TRUE),
		time = c(19, 19.5, 11.2, 11, 7.5, 8, 8.25, 8.4)
		)

	df1.mean.se <- ddply(df1, .(strain), summarise, "sem" = sd(time)/sqrt(length(time)), "mean" = mean(time))
	limits1 <- aes(ymax = mean + sem, ymin=mean - sem)

	f5b <- ggplot(df1.mean.se, aes(strain, mean)) + geom_bar(stat="identity") + geom_errorbar(limits1, width=0.2) + 
	labs(x="Strain", y ="Generation time (h)") + 
	theme_bw(base_size = 8, base_family = "ArialMT") + theme(axis.title.x = element_blank())
}

f5c1 = func5c() {
	df2 <- data.frame(
		strain = c(rep("p2", 3),rep("p20",3), rep("p2",3 ), rep("p20", 3)),
		organ = factor(c(rep("spleen", 6),rep("liver", 6)), levels = c("spleen", "liver"), ordered=TRUE),
		burden = c(50000000000, 50000000000, 50000000000, 125000, 77000, 102000, 97500000000, 10000000000, 50000000000, 660000, 1900000, 1000000)
		)
	df2.mean.se <- ddply(df2, .(strain, organ), summarise, "sem" = sd(burden)/sqrt(length(burden)), "mean" = mean(burden))
	limits2 <- aes(ymax = mean + sem, ymin=mean - sem)	

	f5c1 <- ggplot(df2.mean.se, aes(x=factor(strain), y=mean)) + 
		geom_bar(stat="identity") + geom_errorbar(limits2, width = 0.2) + scale_y_log10() + 
		facet_grid(~organ) + 
		labs(x = "strain", y = "parasite burden")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
		science_theme2 +
		theme(  
	     axis.title.x  = element_blank(),
	     legend.direction = 'horizontal', legend.position = 'top', 
	     strip.background = element_blank(),strip.text.x = element_text(face="bold")) + 
		guides(fill=guide_legend(nrow=1,byrow=FALSE))

}

f5c2 = function() {
	df3 <- data.frame(
		weight = c( 104.3, 125.5, 121.5, 120.8, 109.3, 108.6, 
					121.6, 135.6, 136, 142.3, 130, 128.5, 
					139.2, 147.5, 143.5, 152.8, 129.7, 136.3, 
					127.8, 141.7, 140.2, 160.6, 135.9, 149,
					125.4, 128.8, 126.2, 170, 147.8, 150.9,
					98.4, 118, 107.5, 174.2, 154, 160.7,
					97.6, 97.7, 103.9, 178, 162.4, 162.1
				),
		time = c(rep(0,6),rep(3,6),rep(5,6),rep(7,6),rep(9,6),rep(12,6),rep(13,6)),
		strain = factor(c(rep("p2", 3),rep("p20", 3),rep("p2", 3),rep("p20", 3),rep("p2", 3),rep("p20", 3),rep("p2", 3),rep("p20", 3),rep("p2", 3),rep("p20", 3),rep("p2", 3),rep("p20", 3), rep("p2", 3),rep("p20", 3)), levels=c("p2","p20"), ordered=TRUE)
		)

		df3.mean.se <- ddply(df3, .(time, strain), summarise, "sem" = sd(weight)/sqrt(length(weight)), "mean" = mean(weight))
		limits3 <- aes(ymax = mean + sem, ymin=mean - sem)	

		f5c2 <- ggplot(df3.mean.se, aes(color=strain, y=mean, x=time)) + 
		geom_point() + geom_line(aes(group=strain)) + geom_errorbar(limits3, width=0.2) + 
		labs(x="time (weeks)", y = "Hamster weight (g)") + 
		theme_bw(base_size = 8, base_family = "ArialMT") + theme(legend.title=element_blank(), legend.position = c(0.82,0.6)) +
		guides(color = guide_legend(keywidth = .5, keyheight = .5))

}

f5b = func5b() {


	read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/data/R/C7A9MANXX_1_15nf.filtered.bed") -> df.test.1
	df.test.1$sample = "spleen"
	read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/data/R/C7A9MANXX_4_8nf.filtered.bed") -> df.test.4
	df.test.4$sample = "liver"

	df.liver.spleen = rbind(df.test.1, df.test.4)
	df.liver.spleen$sample = factor(df.liver.spleen$sample)
	colnames(df.liver.spleen) <- c("chromosome", "start", "freq", "sample")
	df.liver.spleen$chromosome = gsub(df.liver.spleen$chromosome, pattern="\\|\\S+", replacement="")
	levels = mixedsort(levels(coverage.melted2$variable)), ordered = TRUE

	df.liver.spleen$chromosome <- factor(array(df.liver.spleen$chromosome), levels = mixedsort(unique(array(df.liver.spleen$chromosome))), ordered = TRUE)
	
	f5d <- ggplot(subset(df.liver.spleen, chromosome %in% c("5", "14", "15", "20", "26")), aes(freq)) +
	geom_density(adjust=1.5) +
	geom_histogram(aes(y=..density..)) +
	#geom_histogram() +
	facet_grid(sample~chromosome) + 
	labs(x = "Frequency", y = "Density", title="") + 
	science_theme2 +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	theme(   axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.y = element_text(face = "bold"),
			legend.position = "none",
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank())
}

postscript("~/plotsave/fig5.eps", width = 6, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	ggdraw() +
	draw_plot(f5a,   0.0773, 	0.6, 	.85, 	.35) +
  	draw_plot(f5b,   0.0755, 	0.28, 	.25,  	.25) +
  	draw_plot(f5c1,  0.34, 	0.28, 	.28,  	.271) + 
  	draw_plot(f5c2,  0.65, 	0.261, 	.317,  	.273) + 
  	draw_plot(f5d,   0.095, 	0, 		0.9,  	.25) + 
	draw_plot_label(c("a", "b", "c"), c(0.04, 0.04, 0.04), c(0.98, 0.58, 0.26), size = 15)
dev.off()

melted.rpkm.t = melt(rpkm[,c(2,8:15)], id.vars=c("CHROMOSOME"))
melted.gpkm.t = melt(gpkm[,c(2,8:15)], id.vars=c("CHROMOSOME"))
medians.melted.t = ddply(melted.rpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})[,1:2]
medians.melted.t$rpkm = ddply(melted.rpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})$V1
medians.melted.t$gpkm = ddply(melted.gpkm.t, .(CHROMOSOME,variable), function(x){median(x$value)})$V1
chrv <- array(medians.melted.t$CHROMOSOME)
medians.melted.t$CHROMOSOME <- factor(chrv, levels=mixedsort(unique(chrv)))
f5a <- ggplot(subset(medians.melted.t, !(CHROMOSOME %in% c("MC", "MI"))), aes(gpkm, rpkm)) + 
geom_point(size = .8) +
geom_abline(slope=1, intercept=-30) + 
geom_point(data=subset(medians.melted.t, (CHROMOSOME %in% c("5", "14", "15", "20", "26", "31"))), aes(color=CHROMOSOME), size = .8) + 
labs(color = "Chromosome", x="Genomic coverage", y="Transcriptomic coverage", title="") + 
theme_bw(base_size = 8, base_family = "ArialMT") + science_theme2 + guides(fill=guide_legend(nrow=3,byrow=TRUE))



################ S11

df.liver.spleen$chrgroup = ""
df.liver.spleen[which(df.liver.spleen$chromosome %in% levels(df.liver.spleen$chromosome)[1:12]),"chrgroup"] <- "first"
df.liver.spleen[which(df.liver.spleen$chromosome %in% levels(df.liver.spleen$chromosome)[13:24]),"chrgroup"] <- "second"
df.liver.spleen[which(df.liver.spleen$chromosome %in% levels(df.liver.spleen$chromosome)[25:36]),"chrgroup"] <- "third"


	fs11_1 <- ggplot(subset(df.liver.spleen, chromosome %in% levels(df.liver.spleen$chromosome)[1:3]), aes(freq)) +
	geom_density(adjust=1.5) +
	geom_histogram(aes(y=..density..)) +
	facet_grid(chrgroup + sample ~ chromosome) + 
	labs(x = "", y = "Density", title="") + 
	science_theme2 +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	theme( plot.margin = unit(c(0,0,0,0), "cm"), panel.margin.y = unit(0.4, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			legend.position = "none")

	fs11_2 <- ggplot(subset(df.liver.spleen, chromosome %in% levels(df.liver.spleen$chromosome)[13:24]), aes(freq)) +
	geom_rect(data = subset(df.liver.spleen, chromosome == "20"),aes(colour = "red", fill = "white"),xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf) +
	geom_density(adjust=1.5) +
	geom_histogram(aes(y=..density..)) +
	#geom_histogram() +
	facet_grid(sample~chromosome) + 
	labs(x = "", y = "Density", title="") + 
	science_theme2 +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	theme( plot.margin = unit(c(0,0,0,0), "cm"), panel.margin.y = unit(0.4, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			legend.position = "none")

	fs11_3 <- ggplot(subset(df.liver.spleen, chromosome %in% levels(df.liver.spleen$chromosome)[25:36]), aes(freq)) +
	geom_density(adjust=1.5) +
	geom_histogram(aes(y=..density..)) +
	#geom_histogram() +
	facet_grid(sample~chromosome) + 
	labs(x = "Frequency", y = "Density", title="") + 
	science_theme2 +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	theme( plot.margin = unit(c(0,0,0,0), "cm"), panel.margin.y = unit(0.4, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			legend.position = "none")

postscript("~/plotsave/figs11.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	ggdraw() +
	draw_plot(fs11_1,   0, .66, 1,  .33) +
  	draw_plot(fs11_2,   0, .33, 1,  .33) +
  	draw_plot(fs11_3,   0, 0, 1,  .33)  	
dev.off()


#

gpkm.passages.m$chrgroup = ""
gpkm.passages.m[which(gpkm.passages.m$Chromosome %in% levels(gpkm.passages.m$Chromosome)[1:18]),"chrgroup"] <- "first"
gpkm.passages.m[which(gpkm.passages.m$Chromosome %in% levels(gpkm.passages.m$Chromosome)[19:36]),"chrgroup"] <- "second"

fs5_1_1 <- ggplot(subset(gpkm.passages.m, chrgroup != "" & variable %in% c("Sp_ama", "p2", "p10", "p20")), aes(Chromosome, value, fill = variable)) + 
geom_boxplot( outlier.shape = NA) + 
coord_cartesian(ylim=c(0,160))  + 
facet_wrap(~chrgroup, scales = "free", nrow = 2) +
scale_fill_manual(values = gg_color_hue(4),labels = c("Splenic amastigote", "p2", "p10", "p20")) +
labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
science_theme2 +
	theme(strip.text.x = element_blank(), strip.background = element_blank(), axis.title.x = element_blank(), legend.direction = 'horizontal', legend.position = 'top') + 
	guides(fill=guide_legend(nrow=1,byrow=FALSE))

#fs5_1_2 <- ggplot(subset(gpkm.passages.m, Chromosome %in% levels(gpkm.passages.m$Chromosome)[19:36] & variable %in% c("Sp_ama", "p2", "p10", "p20")), aes(Chromosome, value, fill = variable)) + 
#geom_boxplot( outlier.shape = NA) + coord_cartesian(ylim=c(0,160))  + 
#scale_fill_manual(values = gg_color_hue(4),labels = c("Splenic amastigote", "P2", "P10", "P20")) +
#labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 8, base_family = "ArialMT") + 
#science_theme2 + theme(legend.position = 'none', plot.margin = unit(c(-0.3,0,-0.3,0), "lines"))

gpkm.clones.m$chrgroup = ""
gpkm.clones.m[which(gpkm.clones.m$Chromosome %in% levels(gpkm.clones.m$Chromosome)[1:18]),"chrgroup"] <- "first"
gpkm.clones.m[which(gpkm.clones.m$Chromosome %in% levels(gpkm.clones.m$Chromosome)[19:36]),"chrgroup"] <- "second"

fs5_3_1 <- ggplot(subset(gpkm.clones.m, chrgroup != "" & variable %in% c("p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")), aes(Chromosome, value, fill = variable)) + 
	geom_boxplot( outlier.shape = NA) + coord_cartesian(ylim=c(0,160)) + 
	facet_wrap(~chrgroup, scales = "free", nrow = 2) +
	labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
	scale_fill_manual(values = gg_color_hue(9),labels = c("p20", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")) +
	theme_bw(base_size = 8, base_family = "ArialMT") + 
	science_theme2 + 
	theme(strip.text.x = element_blank(), strip.background = element_blank(), legend.title  = element_blank(), legend.position = 'bottom', legend.direction = 'horizontal') +
	guides(fill=guide_legend(nrow=1,byrow=FALSE))

#fs5_3_2 <- ggplot(subset(gpkm.clones.m, Chromosome %in% levels(gpkm.clones.m$Chromosome)[19:36] & variable %in% c("p21", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")), aes(Chromosome, value, fill = variable)) + 
#	geom_boxplot( outlier.shape = NA) + coord_cartesian(ylim=c(0,160)) + 
#	labs(x = "Chromosome", y = "Read depth", fill = "Sample", title="")  + theme_bw(base_size = 12, base_family = "Helvetica") + 
#	scale_fill_manual(values = gg_color_hue(9),labels = c("P20", "CL1", "CL3", "CL4", "CL6", "CL7", "CL8", "CL9", "CL10")) +
#	theme_bw(base_size = 8, base_family = "ArialMT") + 
#	science_theme2 + 
#	theme(plot.margin = unit(c(-0.3,0,0,0) "lines"), legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'bottom') + 
#	guides(fill=guide_legend(nrow=1,byrow=FALSE))


postscript("~/plotsave/figs5.eps", width = 8, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	ggdraw() +
	draw_plot(fs5_1_1,	 0,		.63, 	1,  	.40) +
  	draw_plot(fs5_3_1,   0, 	.20, 	1,  	.40) +
  	draw_plot(fs5_2,   	.33, 	0, 	.33,  	.20) +
	draw_plot_label(c("a", "b", "c"), c(0, 0, 0.33), c(1, .6, .2), size = 15)
dev.off()
if(FALSE) {



fs5_2 <- autoplot(coverage.clones.pca,data = coverage.clones, colour = 'sample',shape = TRUE, label = FALSE, label.size = 3) + 
coord_cartesian(xlim=c(-0.4,0.7)) + 
scale_colour_manual(values = gg_color_hue(9)[2:9]) +
labs(x="PC1 (77.92% explained var.)", y ="PC2 (15.30% explained var.)", colour = '') + 
theme_bw(base_size = 8, base_family = "ArialMT")

postscript("~/tmp/newfigs5.eps", width = 8, height = 10, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	ggdraw() +
	draw_plot(fs5_1_1,	 0,		.63, 	1,  	.40) +
  	draw_plot(fs5_3_1,   0, 	.23, 	1,  	.40) +
  	draw_plot(fs5_2,   	.33, 	0, 	.33,  	.20) +
	draw_plot_label(c("a", "b", "c"), c(0, 0, 0.30), c(1, .6, .2), size = 15)
dev.off()
}


coverage.clones <- read.csv("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/data/gene_counts/DNA_09_2014_subclones/coverage.txt", header = TRUE)[3:10,]

row.names(coverage.clones) <- coverage.clones$sample
coverage.clones$sample <- row.names(coverage.clones)

coverage.clones.pca <- prcomp(coverage.clones[1:36], scale. = FALSE)



fs5_2 <- autoplot(coverage.clones.pca,data = coverage.clones, colour = 'samples',shape = TRUE, label = FALSE, label.size = 3) + coord_cartesian(xlim=c(-0.4,0.7)) + labs(x="PC1 (77.92% explained var.)", y ="PC2 (15.30% explained var.)")+ theme_bw(base_size = 8, base_family = "ArialMT")









################################### Figure Supp 3 A+B

df.het.by.somy.temp <- subset(df.het.by.somy, somy %in% c("disomic", "trisomic"))
df.het.by.somy.temp$somy <- "any"
df.het.by.somy.2 <- rbind(df.het.by.somy, df.het.by.somy.temp)


equal_breaks <- function(x){
    q <- round(range(x, na.rm = TRUE)[2] / 100)
    if(q < 3) {
    	breaks <- seq(0, q*100, by = 100)
    } else {
    	breaks <- seq(0,q*100, by = round((q/2))*100)
    }

    names(breaks) <- attr(breaks,"labels")
    breaks
}


fS3A <- ggplot(subset(df.het.by.somy.2, somy %in% c("any", "disomic", "trisomic")), aes(frequency)) + 
		geom_histogram(binwidth=0.01) + 
		facet_grid(chromosome~somy, scales="free_y") + 
		scale_y_continuous(breaks=equal_breaks) +
		scale_x_continuous(breaks = c(0.25,0.5,0.75)) +
		theme(axis.text.x = element_text(angle=90)) +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		science_theme2

		df.somy.counts$chromosome = factor(array(df.somy.counts$chromosome), rev(levels(df.somy.counts$chromosome)), ordered = TRUE)

fS3B <- ggplot(df.somy.counts, aes(chromosome, samples, fill = somy)) + geom_bar(stat = "identity") +
		coord_flip() +
		labs(x= "", y = "Number of samples", fill = "Somy") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		science_theme2 +
		theme(axis.text.y=element_blank(), axis.ticks.y = element_blank())

postscript("~/plotsave/figs3.eps", width = 8, height = 14, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	
	ggdraw() +
	draw_plot(fS3A,	 0,		0, 		0.4,  	1) +
  	draw_plot(fS3B,   0.4, 	0, 	0.6,  	0.985) +
	draw_plot_label(c("a", "b"), c(0, 0.4), c(1, 1), size = 15)

	
dev.off()



######################

read.table("~/projects/kinetoplastids/Ldonovani/data/R/samples_to_snp_count.counts.txt") -> df.snps.field
colnames(df.snps.field) <- c("chromosome", "sites", "length")
vars <- gsub(df.snps.field$chromosome, pattern="\\|\\S+", replacement="")
df.snps.field$chromosome <- factor(vars, levels = unique(mixedsort(vars)), ordered = TRUE)
num_samples <- length(levels(coverage.melted$Isolate))
df.aneuploid.field <- ddply(coverage.melted, c("variable"), function(x){ length(which(x$value > 1.2)) / num_samples })[1:36,]
chrreplace <- gsub(df.aneuploid.field$variable, pattern="Ld", replacement="")
df.aneuploid.field$variable <- factor(chrreplace, levels = , ordered = TRUE)

df.snps.field$anefreq <- rep(0,38)
df.snps.field[3:38,"anefreq"] <- df.aneuploid.field$V1



postscript("~/plotsave/fig6.eps", width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")

ggplot(df.snps.field, aes(length/1000000, sites, label=chromosome)) + geom_point(size = 6, shape = 16, aes(colour = anefreq)) + geom_text(size=3) + 
geom_abline(color="blue", slope = 0.0004*1000000) + geom_abline(color="green", slope = 0.0007*1000000) + geom_abline(color="red", slope = 0.0018*1000000) +
scale_x_continuous(breaks=c(1,2),labels=c("1","2")) +
#scale_shape_discrete(solid=T) +
scale_colour_gradient(low = "white", high = "#EE3B3B") +
labs(colour = "Aneuploidy\nfrequency", x="Chromosome length (Mb)", y="# Het variant sites") + 
theme_bw(base_size = 8, base_family = "ArialMT") + science_theme2
dev.off()

#####################


df.field.piles <- read.table("~/projects/kinetoplastids/Ldonovani/data/all.pile.Rpile.depth.1.6.smiley.new.txt")
colnames(df.field.piles) <- c("chromosome", "start", "frequency", "sample", "depth")
newchr <- gsub(array(df.field.piles$chromosome), pattern="Ld", replacement = "")
df.field.piles$chromosome <- factor(newchr, levels = mixedsort(unique(newchr)), ordered = TRUE)

select.chromosomes1 <- c("5", "12")
select.samples1 <- row.names(coverage)[which(coverage$Ld12 > 1.15 & coverage$Ld5 > 1.15)]
fS4A <- ggplot(subset(df.field.piles, chromosome %in% select.chromosomes1 & sample %in% select.samples1), aes(frequency)) + geom_histogram() + 
facet_grid(chromosome~sample) + 
theme_bw(base_size = 8, base_family = "ArialMT") +
theme(axis.text.x = element_text(angle = 90), panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"))

select.chromosomes2 <- c("6", "12")
select.samples2 <- row.names(coverage)[which(coverage$Ld12 > 1.15 & coverage$Ld6 > 1.15)]
fS4B <- ggplot(subset(df.field.piles, chromosome %in% select.chromosomes2 & sample %in% select.samples2), aes(frequency)) + geom_histogram() + 
facet_grid(chromosome~sample) +
theme_bw(base_size = 8, base_family = "ArialMT") +
theme(axis.text.x = element_text(angle = 90), panel.grid.major = element_line(size = 0.5, color = "grey"), 
    axis.line = element_line(size = 0.7, color = "black"))

postscript("~/plotsave/figs4.eps", width = 6, height = 4, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	
	ggdraw() +
	draw_plot(fS4A,	 0,		0.5, 		1,  	0.5) +
  	draw_plot(fS4B,   0, 	0, 	1,  	0.5) +
	draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.5), size = 15)

	
dev.off()

################# Genome coverage bins ############################

df.field.chrbins <- read.table("/users/cn/pprieto/projects/kinetoplastids/Ldonovani/analysis/2015-12-Hideo/mapping_and_calling/mapping_data/last_bpk/merged.chrbin.depth.noempty.txt")
rownames(df.field.chrbins) <- array(df.field.chrbins$V2)
data.prop.1 <- df.field.chrbins[, 3:382]
df.field.chrbins.m <- melt(df.field.chrbins, id.vars = c("V1", "V2"))
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
chrgroups = {}
for(c in 1:38){ chrgroups = c(chrgroups, rep(c,10*length(row.names(df.field.chrbins))))}
df.field.chrbins.m$chrgroup = chrgroups

postscript("~/plotsave/figs1.eps", width = 14, height = 32, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")

ggplot(subset(df.field.chrbins.m, !(chrgroup %in% c(37:38))), aes(x=variable, y=V2, fill=value*2)) + geom_tile() + 
facet_grid(.~chrgroup, scales = "free") +
scale_fill_gradientn(limits = c(1.6,4.5),  colours = c("yellow", "green", "blue")) +
labs(x = "Chromosome", y = "Sample", fill = "Depth") +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

############## Correlation plot ##########################


# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

reorder_cormat <- function(cormat){
	# Use correlation between variables as distance
	dd <- as.dist((1-cormat)/2)
	hc <- hclust(dd)
	cormat <-cormat[hc$order, hc$order]
}

#cormat <- round(cor(coverage_reformatted[1:37]),2)
cmat <- round(cor(coverage),2)

cmat.r <- reorder_cormat(cmat)
upper_tri <- get_upper_tri(cmat.r)

melted_cormat <- melt(as.matrix(upper_tri), na.rm = TRUE)

melted_cormat$X1 <- factor(array(melted_cormat$X1), levels = dimnames(cmat.r)[[1]], ordered = TRUE)
melted_cormat$X2 <-  factor(array(melted_cormat$X2), levels = dimnames(cmat.r)[[1]], ordered = TRUE)

#melted_cormat <- reshape2::melt(cormat)[reshape2::melt(upper.tri(cormat))[, 3], ]

ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value))+
 geom_tile(color = "white") +
  scale_fill_gradient2(na.value = 'white', low = "gold", high = "red4", mid = "orange", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 8, hjust = 1))+
 coord_fixed()

cor.mtest <- function(mat, conf.level = 0.95){
  mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    diag(lowCI.mat) <- diag(uppCI.mat) <- 1
    for(i in 1:(n-1)){
        for(j in (i+1):n){
            tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
            p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
            lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
            uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
        }
    }
    return(list(p.mat, lowCI.mat, uppCI.mat))
}

colnames(coverage) <- gsub(colnames(coverage), pattern = "Ld", replacement = "")
coverage.sub <- round(coverage[rowSums(is.na(coverage)) == 0,1:36], digits=1)
res1 <- cor.mtest(coverage.sub,0.95)
M <- cor(coverage.sub)
diag(M) <- NA
## specialized the insignificant value according to the significant level
#corrplot.mixed(M, na.label = "o", insig = "blank", tl.col = "black", order="hclust")

corrplot(M, p.mat = res1[[1]], sig.level=0.01, na.label = "o", insig = "blank", tl.col = "black", order="hclust", tl.cex = 1, cl.cex = 1)

library(gridGraphics)
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

g <- grab_grob()

plot <- arrangeGrob(g, ncol=1)


postscript("~/plotsave/figs2new.eps", width = 16, height = 16, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
corrplot(M, p.mat = res1[[1]], insig = "blank")
dev.off()
ggheatmap +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))

  dev.off()
  

  
  
postscript("~/tmp/new_fig1.eps", width = 8, height = 6, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
	
png("~/tmp/newfig1.png", width = 800, height = 800)
ggdraw() +
  	draw_plot(plot, .3, .52, .5, .5) +
	draw_plot(new_f1a, 0, 0, 0.2, 1) +
	draw_plot(new_f1c.1, .25, 0.35, .32, 0.2) +
  	draw_plot(new_f1c.2, (0.25 + 0.3), 0.35, 0.2, 0.18) +
  	draw_plot(new_f1c.3, (0.25 + 0.5), 0.35, 0.22, 0.18) +
  	draw_plot(new_f1d, .25, 0, .7, .35) +
  	draw_plot_label(c("a", "b", "c", "d"), c(0, .25, .25, .25), c(1, 1, .55, .35), size = 15)
dev.off()