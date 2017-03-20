# To obtain the dataset:
cat all/*_real_dominant_disomic_nsamples_results.txt.counts | awk '{OFS="\t"; print $1"_"$2,$4,$5,$7,$8,$9,$10}' | sort -k 1,1 > all/real_dominant_disomic_nsamples_results.txt.counts
cat all/*_real_dominant_trisomic_nsamples_results.txt.counts.2 | awk '{OFS="\t"; print $1"_"$2,$4,$5,$7,$8,$9,$10}' | sort -k 1,1 > all/real_dominant_trisomic_nsamples_results.txt.counts

join all/real_dominant_disomic_nsamples_results.txt.counts all/real_dominant_trisomic_nsamples_results.txt.counts > real_dominant_somic_nsamples_results.txt.counts
# Then change _ for white space and rename BS chars in chromosme names

# Libraries
library(foreach)
library(doMC)
registerDoMC(cores = 4)

# Version to build distributinos of log fold changes
read.table("../contrast_alleles/real_dominant_somic_nsamples_results.txt.counts") -> df.dominante
colnames(df.dominante) <- c("chr", "start", "dicount", "ditotal" , "da1", "da2", "disumfreq1", "disumfreq2", "tricount", "tritotal", "ta1", "ta2", "trisumfreq1", "trisumfreq2")
df.dominante$diffmaxfreq = 0
for(i in 1:dim(df.dominante)[1]) {
	maxal = 0
	maxdifreq = 0
	maxtrifreq = 0

	if(df.dominante[i,'trisumfreq1'] > df.dominante[i,'trisumfreq1']) {
		maxal = df.dominante[i,'ta1']
		maxtrifreq = df.dominante[i,'trisumfreq1']
	} else {
		maxal = df.dominante[i,'ta2']
		maxtrifreq = df.dominante[i,'trisumfreq2']
	}

	if(df.dominante[i,'da1'] == maxal) {
		maxdifreq = df.dominante[i,'disumfreq1']
	} else if(df.dominante[i,'da2'] == maxal) {
		maxdifreq = df.dominante[i,'disumfreq2']
	}
	df.dominante[i,'diffmaxfreq'] <- (log2(maxtrifreq) / log2(maxdifreq))
}

png("~/tmp/log2foldc_2.png", 300, 2000, pointsize=20)
ggplot(subset(df.dominante, ditotal > 2 & tritotal > 2), aes(diffmaxfreq)) + 
geom_histogram(aes(y=..density..)) + geom_density(fill=NA, colour="black") +
facet_grid(chr~., scales = "free") +
labs(x = "Log fold change")
dev.off()

## New table
read.table("real_dominant_somic_nsamples_results.txt.counts") -> df.dominante
colnames(df.dominante) <- c("chr", "start", "dicount", "ditotal" , "disumfreq1", "disumfreq2", "tricount", "tritotal", "trisumfreq1", "trisumfreq2")

df.dominante.old <- df.dominante
df.dominante <- subset(df.dominante, ditotal > 3 & tritotal > 3 )
df.melt1 <- melt(df.dominante, id.vars = c("chr", "start", "dicount", "ditotal" , "tricount", "tritotal"), measure.vars = c("disumfreq1", "disumfreq2"))
df.melt1$variable <- "disomic"
df.melt2 <- melt(df.dominante, id.vars = c("chr", "start", "dicount", "ditotal" , "tricount", "tritotal"), measure.vars = c("trisumfreq1", "trisumfreq2"))
df.melt2$variable <- "trisomic"

df.sumfreq.melt <- rbind(df.melt1, df.melt2)
mixedsort(unique(array(df.dominante$chr))) -> chrs
factor(df.dominante$chr, levels = chrs, ordered = TRUE) -> df.dominante$chr

png("~/tmp/dominant_sumprofiles_filter.png", 600, 1000, pointsize=20)
ggplot(df.sumfreq.melt, aes(value)) + geom_density() + facet_grid(chr ~ variable, scales = "free")
dev.off()

png("~/tmp/dominant_sumhist_filter_onlydisomyfilter.png", 600, 1000, pointsize=20)
ggplot(df.sumfreq.melt, aes(value)) + geom_histogram() + facet_grid(chr ~ variable, scales = "free")
dev.off()

# Previous stuff

read.table("real_dominant_filtered_somic_nsamples_results.txt") -> df.dominante
colnames(df.dominante) <- c("chr", "start", "dicount", "ditotal", "ditype" , "tricount", "tritotal", "tritype")

df.dominante$trifreq <- df.dominante$tricount / df.dominante$tritotal
df.dominante$difreq <- df.dominante$dicount / df.dominante$ditotal

df.dominante$freqdiff <- df.dominante$trifreq - df.dominante$difreq

mixedsort(unique(array(df.dominante$chr))) -> chrs
factor(df.dominante$chr, levels = chrs, ordered = TRUE) -> df.dominante$chr

melt(df.dominante, measure.vars=c("difreq", "trifreq")) -> df.dominante.m


ggplot(df.dominane, aes(difreq, trifreq, color = chr)) + geom_histogram() 

ggplot(df.dominante, aes(difreq, trifreq, color = chr)) + geom_histogram() + facet_grid(chr~.)

ggplot(df.dominante, aes(trifreq-difreq, fill = chr)) + geom_density(alpha = .4) 

png("~/tmp/dominant_profiles_filt.png", 600, 1000, pointsize=20)
ggplot(df.dominante.m, aes(value))+ geom_histogram() + facet_grid(chr~variable, scales = "free")
dev.off()
ggplot(df.dominante.m aes(tritotal))+ geom_density() + facet_grid(chr~.)

df.dominante$ditotalcut <-  cut(df.dominante$ditotal, breaks = c(0,3,204))


# Subset the sites that seem to be supported by at least more than X disomic samples...
# Compared against the unique sites in disomics
png("~/tmp/dominant_diff_cuts.png", 600, 1000, pointsize=20)
ggplot(df.dominante, aes(trifreq-difreq, fill = chr)) + geom_density(alpha = .4) + facet_grid(chr~ditotalcut, scales = "free")
dev.off()

png("~/tmp/dominant_diff_hist_cuts.png", 600, 1000, pointsize=20)
ggplot(df.dominante, aes(trifreq-difreq, fill = chr)) + geom_histogram(alpha = .4) + facet_grid(chr~ditotalcut, scales = "free")
dev.off()


df.dominante.subdifferent <- subset(df.dominante, !(trifreq == 1 & difreq == 1))

png("~/tmp/dominant_subdiff_hist_cuts.png", 600, 1000, pointsize=20)
ggplot(subset(df.dominante.subdifferent), aes(trifreq-difreq, fill = chr)) + geom_histogram(alpha = .4) + facet_grid(chr~., scales = "free")
dev.off()


###### Tomachina

	df.het <- read.table("all.pile.Rpile.depth.1.6.smiley.allele.new.txt",header=F)
	colnames(df.het) <- c("chromosome", "start","frequency", "alelle","name", "depth")
	newchrs <- gsub(pattern="Ld", df.het$chromosome, replacement = "")
	df.het$chromosome = factor( newchrs, levels= mixedsort(unique(newchrs)), ordered=TRUE )

	df.het.by.somy <- ddply(df.het, .(chromosome), function(d) {
		chr = unique(array(d$chromosome))
		if(chr %in% c("MC", "MI")) {
			return()
		}
		array_coverage 		= get(paste("Ld", chr, sep = ""), coverage)
		array_isolates 		= array(rownames(coverage))

		monosomic_samples   = array_isolates[array_coverage < 0.85 & !is.na(array_coverage)]
		disomic_samples     = array_isolates[array_coverage > 0.85 & array_coverage < 1.15 & !is.na(array_coverage)]
		trisomic_samples 	= array_isolates[array_coverage > 1.3 & array_coverage < 1.7 & !is.na(array_coverage)]
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


# 
cmp = {}
for(k in 1:dim(df.dominante)[1]) { 



benchmark(
    foreach(n = 1:100, .combine = c) %do% myfunction(n),
    
)


plot(x, y, xlab = "Aggregate ratio frequency", ylab = "Individual ratio frequency", main = "LALA", type = "n")
points(x, y, col = "black", pch = 16)

foreach(n = 1:10, .combine = c) %dopar% myfunction(n)
	
myfunction <- function(k) {
	maxal = 5
	if(df.dominante[k,"trisumfreq1"] > df.dominante[k,"trisumfreq2"]) {
			maxal = df.dominante[k,"ta1"]
	} else {
			maxal = df.dominante[k,"ta2"]
	}
	disomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "disomic" & alelle == maxal)
	trisomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "trisomic" & alelle == maxal)
	
	if( (dim(disomic)[1] != 0) && (dim(trisomic)[1]) != 0) {
		acc <- rep(NA, (dim(disomic)[1] * dim(trisomic)[1])*3)
		n = 1
		for(i in 1:dim(disomic)[1]) {
			for(j in 1:dim(trisomic)[1]) {
				acc[n] <- df.dominante[k,"diffmaxfreq"]
				acc[n + 1] <- (log2(trisomic[j,"frequency"]) / log2(disomic[i,"frequency"]))
				acc[n + 2] <- df.dominante[k,"chr"]
				n <- n + 3
			}
		}
		acc
	}
}

myavgfunction <- function(k) {
	maxal = 5
	if(df.dominante[k,"trisumfreq1"] > df.dominante[k,"trisumfreq2"]) {
			maxal = df.dominante[k,"ta1"]
	} else {
			maxal = df.dominante[k,"ta2"]
	}
	disomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "disomic" & alelle == maxal)
	trisomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "trisomic" & alelle == maxal)
	
	if( (dim(disomic)[1] != 0) && (dim(trisomic)[1]) != 0) {
		acc <- rep(NA, (dim(trisomic)[1])*3)
		n = 1
		for(j in 1:dim(trisomic)[1]) {
			acc[n] <- df.dominante[k,"diffmaxfreq"]
			acc[n + 1] <- (log2(trisomic[j,"frequency"]) / log2(mean(disomic$frequency, na.rm = TRUE)))
			acc[n + 2] <- df.dominante[k,"chr"]
			n <- n + 3
		}
		acc
	}
}

myCedfunction <- function(k) {
	maxal = 5
	disumfreq = 0
	if(df.dominante[k,"trisumfreq1"] > df.dominante[k,"trisumfreq2"]) {
			maxal = df.dominante[k,"ta1"]
	} else {
			maxal = df.dominante[k,"ta2"]
	}

	if(df.dominante[k,"da1"] == maxal) {
		disumfreq = df.dominante[k,"disumfreq1"]
	}
	if(df.dominante[k,"da2"] == maxal) {
		disumfreq = df.dominante[k,"disumfreq2"]
	}
	
	trisomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "trisomic" & alelle == maxal)
	
	if( (dim(trisomic)[1]) != 0) {
		acc <- rep(NA, (dim(trisomic)[1])*2)
		n = 1
		for(j in 1:dim(trisomic)[1]) {
			acc[n] <- (log2(trisomic[j,"frequency"]) / log2(disumfreq))
			acc[n + 1] <- df.dominante[k,"chr"]
			n <- n + 2
		}
		acc
	}
}

myPabfunction <- function(k) {
	maxal = 5
	disumfreq = 0
	trisumfreq = 0
	if(df.dominante[k,"trisumfreq1"] > df.dominante[k,"trisumfreq2"]) {
			maxal = df.dominante[k,"ta1"]
			trisumfreq = df.dominante[k,"trisumfreq1"]
	} else {
			maxal = df.dominante[k,"ta2"]
			trisumfreq = df.dominante[k,"trisumfreq2"]
	}

	if(df.dominante[k,"da1"] == maxal) {
		disumfreq = df.dominante[k,"disumfreq1"]
	}
	if(df.dominante[k,"da2"] == maxal) {
		disumfreq = df.dominante[k,"disumfreq2"]
	}
	
	trisomic = subset(df.het.by.somy, chromosome == df.dominante[k,"chr"] & start == df.dominante[k,"start"] & somy == "trisomic" & alelle == maxal)
	
	if( (dim(trisomic)[1]) != 0) {
		acc <- rep(NA, (dim(trisomic)[1])*3)
		n = 1
		for(j in 1:dim(trisomic)[1]) {
			#acc[n] <- (log2(trisomic[j,"frequency"]) / log2(disumfreq))
			acc[n] <- (log2(trisomic[j,"frequency"]/ disumfreq))
			acc[n + 1] <- (log2(trisumfreq/disumfreq))
			#acc[n + 1] <- (log2(trisumfreq /disumfreq))
			acc[n + 2] <- df.dominante[k,"chr"]
			n <- n + 3
		}
		acc
	}
}

plotarr <- foreach(n = 1:dim(df.dominante)[1], .combine = c) %dopar% myfunction(n)
df.res <- data.frame(aggregated = plotarr[seq(from=1,to=length(plotarr), by = 3)], individual = plotarr[seq(from=2,to=length(plotarr), by = 3)], chr = plotarr[seq(from=3,to=length(plotarr), by = 3)])
df.res$chr <- factor(df.res$chr, levels = mixedsort(unique(df.dominante$chr)), ordered = TRUE)

png("~/tmp/aggregated_vs_individual.png", 400, 4000, pointsize=20)
ggplot(df.res, aes(aggregated, individual)) + geom_point() + facet_grid(chr~.) + labs(y = "Individual frequency fold change", x = "Aggregated frequency fold change") + theme_grey()
dev.off()

plotarravg <- foreach(n = 1:dim(df.dominante)[1], .combine = c) %dopar% myavgfunction(n)
df.res.avg <- data.frame(aggregated = plotarravg[seq(from=1,to=length(plotarravg), by = 3)], individual = plotarravg[seq(from=2,to=length(plotarravg), by = 3)], chr = plotarravg[seq(from=3,to=length(plotarravg), by = 3)])
df.res.avg$chr <- factor(df.res.avg$chr, levels = mixedsort(unique(df.dominante$chr)), ordered = TRUE)

png("~/tmp/aggregated_vs_individual_avgdisomic_ver.png", 400, 8000, pointsize=20)
ggplot(df.res.avg, aes(aggregated, individual)) + geom_point() + facet_grid(chr~.) + labs(y = "Individual frequency fold change", x = "Aggregated frequency fold change") + theme_grey()
dev.off()

plotarrced <- foreach(n = 1:dim(df.dominante)[1], .combine = c) %dopar% myCedfunction(n)
df.res.ced <- data.frame(individual = plotarrced[seq(from=1,to=length(plotarrced), by = 2)], chr = plotarrced[seq(from=2,to=length(plotarrced), by = 2)])
df.res.ced$chr <- factor(df.res.ced$chr, levels = mixedsort(unique(df.dominante$chr)), ordered = TRUE)

png("~/tmp/aggregated_vs_individual_sumdisomic_ver.png", 400, 8000, pointsize=20)
ggplot(df.res.ced, aes(individual)) + geom_density() + facet_grid(chr~.) + labs(x = "Individual frequency fold change") + theme_grey()
dev.off()

plotarrpab <- foreach(n = 1:dim(df.dominante)[1], .combine = c) %dopar% myPabfunction(n)
df.res.pab <- data.frame(individual = plotarrpab[seq(from=1,to=length(plotarrpab), by = 3)], aggregated = plotarrpab[seq(from=2,to=length(plotarrpab), by = 3)], chr = plotarrpab[seq(from=3,to=length(plotarrpab), by = 3)])
df.res.pab$chr <- factor(df.res.pab$chr, levels = mixedsort(unique(df.dominante$chr)), ordered = TRUE)

png("~/tmp/aggregated_vs_individual_aggregated_ver.png", 400, 8000, pointsize=20)
ggplot(df.res.pab, aes(aggregated, individual)) + geom_point() + facet_grid(chr~.) + labs(y = "Individual frequency fold change", x = "Aggregated frequency fold change") + theme_grey()
dev.off()

png("~/tmp/aggregated_dist_0.png", 300, 4000, pointsize=20)
ggplot(df.res.pab, aes(aggregated)) + geom_density() + facet_grid(chr~.)
dev.off()

png("~/tmp/individual_dist_0.png", 300, 4000, pointsize=20)
ggplot(df.res.pab, aes(individual)) + geom_density() + facet_grid(chr~.)
dev.off()

df.res.pab.sum <- ddply(df.res.pab, .(chr), function(x) { c(mean(x$individual), mean(x$aggregated))} )
colnames(df.res.pab.sum) <- c("chr", "individual", "aggregated")

lm.sum <- summary(lm(individual~aggregated,data=df.res.pab.sum))
r <- lm.sum$r.squared
f <- lm.sum$fstatistic
p <- pf(f[1],f[2],f[3],lower.tail=F)

df.aneuploid.field <- ddply(coverage.melted, c("variable"), function(x){ length(which(x$value > 1.2 & x$value < 1.8)) / num_samples })[1:36,]
chrreplace <- gsub(df.aneuploid.field$variable, pattern="Ld", replacement="")
df.aneuploid.field$variable <- factor(chrreplace, levels = , ordered = TRUE)
df.res.pab.sum$anefreq <- df.aneuploid.field$V1

df.res.pab.sum.ane <- merge(df.res.pab.sum, df.aneuploid.field, by.x = "chr", by.y = "variable")
colnames(df.res.pab.sum.ane) <- c("chr", "individual", "aggregated", "anefreq")

png("~/tmp/aggregated_vs_individual_means_ver.png", 800, 800, pointsize=20)
new_f1d <-  ggplot(df.res.pab.sum.ane, aes(individual, aggregated, label = chr, fill = anefreq)) + 
geom_label() + 
geom_smooth(method = "lm", se = TRUE, colour = "blue") +
annotate("text", label = paste("P-value =", format(signif(p), scientific = TRUE), sep = " "), x = 0.2, y = 0.4, size = 4, colour = "blue") +
annotate("text", label = paste("R^2 ==", signif(r, digits = 4), sep = " "), x = 0.2, y = 0.38, size = 4, colour = "blue", parse=TRUE) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + 
labs(x="Dispersion", y ="Segregation", fill = "Trisomy frequency" ) +
theme_classic(base_size = 12, base_family = "ArialMT")
dev.off()


### VXZ ###
# sorting chromosomes
# For ratio distributions
ordered_levels <- c(array(df.res.pab.sum[order(df.res.pab.sum$aggregated, decreasing = TRUE), "chr"]), c(1:36)[!(c(1:36) %in% df.res.pab.sum$chr)])
df.res.pab$chr <- factor(df.res.pab$chr, levels = ordered_levels, ordered = TRUE)
# For aggregated frequencies
df.merged.frequencies.piles$chromosome <- factor(df.merged.frequencies.piles$chromosome, levels = ordered_levels, ordered = TRUE)

subset.df.merged.frequencies.piles <- subset(df.merged.frequencies.piles, somy %in% c("disomic", "trisomic"))
subset.df.merged.frequencies.piles$somy <- factor(array(subset.df.merged.frequencies.piles$somy), levels = c("disomic", "trisomic"))

merged_p <- ggplot(subset.df.merged.frequencies.piles, aes(V1)) + 
	geom_histogram() + 
	facet_grid(chromosome~somy, scales = "free", drop = FALSE) + 
	labs(x = "frequency") +
	theme(strip.background = element_blank(), strip.text.x = element_text(face="bold") ,strip.text.y = element_blank(),legend.position = 'none')


aggr_p <- ggplot(df.res.pab) + stat_density(aes(x=aggregated, y=..scaled..), geom="line") + facet_grid(chr~., drop = FALSE) + labs(x="dispersion", y = "density") + geom_vline(aes(xintercept = mean(subset(df.res.pab, chr == 29)$aggregated))) + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = 'none')

ind_p <- ggplot(df.res.pab) + stat_density(aes(x=individual, y=..scaled..), geom="line") + facet_grid(chr~., drop = FALSE) + labs(x = "segregation", y = "density") + geom_vline(aes(xintercept = mean(subset(df.res.pab, chr == 29)$individual))) + theme(strip.background = element_blank(),legend.position = 'none', strip.text.y = element_text(face="bold"))

png("~/tmp/newfigs1_0.png", width = 600, height = 1800)
	ggdraw() +
	draw_plot(merged_p, 0, 0, (1/4)*2, 1) +
  	draw_plot(aggr_p, (0 + (1/4)*2)-0.03, 0.002, (1/4), 0.985) +
  	draw_plot(ind_p, (0 + (1/4)*3)-0.05, 0.002, (1/4)+0.05, 0.985) +
  	draw_plot_label(c("a","b", "c"), c(0,(((1/4)*2)),((1/4)*3)), c(1,1,1), size = 15)
dev.off()

subset.df.merged.frequencies.piles.1c <- subset(df.merged.frequencies.piles, somy %in% c("disomic", "trisomic") & chromosome %in% c(33,35))

new_f1c.1 <- ggplot(subset.df.merged.frequencies.piles.1c, aes(V1)) + 
	geom_histogram() + 
	facet_grid(chromosome~somy, scales = "free") + 
	labs(x = "frequency") +
	theme(strip.background = element_blank(), strip.text.x = element_text(face="bold") ,strip.text.y = element_blank(),legend.position = 'none')
new_f1c.2 <- ggplot(subset(df.res.pab, chr %in% c("35", "5"))) +  stat_density(aes(x=aggregated, y=..scaled..), geom="line") + facet_grid(chr~.) + geom_vline(aes(xintercept = mean(subset(df.res.pab, chr == 29)$aggregated))) + labs(x="dispersion", y = "density") + theme(strip.background = element_blank(), strip.text.y = element_blank(),legend.position = 'none')
new_f1c.3 <- ggplot(subset(df.res.pab, chr %in% c("35", "5"))) +  stat_density(aes(x=individual, y=..scaled..), geom="line") + facet_grid(chr~.) + geom_vline(aes(xintercept = mean(subset(df.res.pab, chr == 29)$individual))) + scale_x_continuous(breaks=seq(-1,1, by = 0.6)) + labs(x = "segregation", y = "density") + theme(strip.background = element_blank(),legend.position = 'none', strip.text.y = element_text(face="bold"))
