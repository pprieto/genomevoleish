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
library(moments)
library(cluster)
library(lfda)
library(ggpmisc)
library(ggrepel)

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

#########

	df.het <- read.table("all.pile.Rpile.depth.1.6.smiley.new.txt",header=F)
	colnames(df.het) <- c("chromosome", "start","frequency","name", "depth")
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
	# Another option to filter by number of samples
	df.allele.std.0 <- read.table("../contrast_alleles/all/all_real_nsamples_results.txt")
	colnames(df.allele.std.0) <- c("chromosome", "start", "allele", "frequency", "samples")
	df.allele.std <- subset(df.allele.std.0, samples > 1)
	df.allele.std$somy <- "std"
	df.het.by.somy.plot <- rbind(df.het.by.somy[,c("chromosome", "start", "frequency", "somy")], df.allele.std[,c("chromosome", "start", "frequency", "somy")])

	

############ To select 
	
read.table("samples_to_snp_count.counts.txt") -> df.snps.field
colnames(df.snps.field) <- c("chromosome", "sites", "length")
vars <- gsub(df.snps.field$chromosome, pattern="\\|\\S+", replacement="")
df.snps.field$chromosome <- factor(vars, levels = unique(mixedsort(vars)), ordered = TRUE)
num_samples <- length(levels(coverage.melted$Isolate))
df.aneuploid.field <- ddply(coverage.melted, c("variable"), function(x){ length(which(x$value > 1.2 & x$value < 1.8)) / num_samples })[1:36,]
chrreplace <- gsub(df.aneuploid.field$variable, pattern="Ld", replacement="")
df.aneuploid.field$variable <- factor(chrreplace, levels = , ordered = TRUE)

df.snps.field$anefreq <- rep(0,38)
df.snps.field[3:38,"anefreq"] <- df.aneuploid.field$V1

df.chr_trisomic = ddply(subset(df.somy.counts, (somy == "trisomic") & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)})
colnames(df.chr_trisomic) <- c("chromosome", "trisomic")

df.chr_sites4std = ddply(df.allele.std, .(chromosome), function(d){dim(d)[1]})
colnames(df.chr_sites4std) <- c("chromosome", "sites4std")

df.chr_surfaceratio = ddply(subset(df.het.by.somy.plot, somy == "trisomic" ), .(chromosome), function(d) {
		den = density(d$frequency)
		tri_area =  sum( den$y[ which( ( den$x > 0.2 & den$x < 0.4 ) | ( den$x > 0.6 & den$x < 0.8 ) ) ] )
		di_area = sum( den$y[ which( (den$x > 0.4 & den$x < 0.6) ) ] )
		tri_area / di_area
	}
)
colnames(df.chr_surfaceratio) <- c("chromosome", "surface_ratio")

df.chr_diffpeaks = ddply( subset(df.het.by.somy.plot, somy == "trisomic"), .(chromosome), function(d){
	den = density(d$frequency, n = 300)
	peaks <- {}
	peakvals <- {}
	for(i in 2:(length(den$y)-1)) {
		k = i - 1
		j = i + 1

		while( (den$y[k] == den$y[i]) && k > 1) {
			k = k - 1
		}

		while( (den$y[j] == den$y[i]) && (j < length(den$y)) ) {
			j = j + 1
		}
		if( den$y[k] < den$y[i] && den$y[i] > den$y[j]) {
			peaks <- c(peaks, i)
			peakvals <- c(peakvals, den$y[i])
		}
	}
	if(length(peaks) == 0 ) {
		return(0)
	}
	if(length(peaks) == 2) {
		return(den$x[peaks[2]] - den$x[peaks[1]])
	}
	maxpeaks = which(round(peakvals, 5) == round(max(peakvals), 5))
	if(length(maxpeaks) == 2) {
		return(den$x[peaks[maxpeaks[2]]] - den$x[peaks[maxpeaks[1]]])
	} else {
		return(0)
	}
})
colnames(df.chr_diffpeaks) <- c("chromosome" , "diff_peaks")

df.meanstd = ddply(subset(df.het.by.somy.plot, somy == "std" ), .(chromosome),  function(d) {
	mean(d$frequency, na.rm = TRUE)
	}
)
colnames(df.meanstd) <- c("chromosome" , "mean_std")

df.medstd = ddply(subset(df.het.by.somy.plot, somy == "std" ), .(chromosome),  function(d) {
	median(d$frequency, na.rm = TRUE)
	}
)
colnames(df.medstd) <- c("chromosome" , "median_std")

df.stdskewness = ddply(subset(df.het.by.somy.plot, somy == "std" ), .(chromosome),  function(d) {
	skewness(d$frequency, na.rm = TRUE)
	}
)
colnames(df.stdskewness) <- c("chromosome" , "std_skewness")

df.stdkurtosis = ddply(subset(df.het.by.somy.plot, somy == "std" ), .(chromosome),  function(d) {
	kurtosis(d$frequency, na.rm = TRUE)
	}
)
colnames(df.stdkurtosis) <- c("chromosome" , "std_kurtosis")

colnames(df.aneuploid.field) <- c("chromosome" , "anefreq" )

df.merged <- Reduce(function(x,y) {merge(x,y, by = "chromosome")}, list(df.chr_trisomic, df.chr_sites4std,df.chr_surfaceratio, df.chr_diffpeaks, df.meanstd, df.medstd, df.stdskewness, df.stdkurtosis, df.aneuploid.field))

# Median case
#df.chr_diffpeaks_surfaceratio = merge( df.chr_diffpeaks , df.chr_surfaceratio, by = "chromosome" )
#colnames(df.chr_diffpeaks) <- c("chromosome" , "diffpeaks")

#df.chr_surfaceratio_medstd = merge( df.chr_surfaceratio, df.medstd, by = "chromosome" )
#colnames(df.chr_surfaceratio_medstd) <- c("chromosome" , "surface_ratio", "median_std")

# Skewness case
#df.chr_surfaceratio_stdskewness = merge( df.chr_surfaceratio, df.stdskewness, by = "chromosome" )
#colnames(df.chr_surfaceratio_stdskewness) <- c("chromosome" , "surface_ratio", "std_skewness")

# Kurtosis case
#df.chr_surfaceratio_stdkurtosis = merge( df.chr_surfaceratio, df.stdkurtosis, by = "chromosome" )
#colnames(df.chr_surfaceratio_stdkurtosis) <- c("chromosome" , "surface_ratio", "std_kurtosis")

#df.chr_surfaceratio_measures = merge( df.chr_surfaceratio_medstd, merge( df.chr_surfaceratio_stdkurtosis, df.chr_surfaceratio_stdskewness, by = c("chromosome", "surface_ratio")), by = c("chromosome", "surface_ratio"))
#df.chr_surfaceratio_stdmeasures_anefreq = merge( df.chr_surfaceratio_measures, df.aneuploid.field, by = "chromosome")


df.chr_surfaceratio_stdmeasures_anefreq.m = melt(df.merged, id.vars = c("chromosome", "diff_peaks", "surface_ratio", "anefreq"), measure.vars= c("mean_std", "median_std", "std_kurtosis", "std_skewness"))
#df.chr_surfaceratio_stdmeasures_anefreq.m = melt(df.chr_surfaceratio_stdmeasures_anefreq, id.vars = c("chromosome", "surface_ratio", "anefreq"), measure.vars= c("median_std", "std_kurtosis", "std_skewness"))



# Still exploring here
# No need to remove the tetrasomic samples!!!!
# Chromosomes to avoid due to lack of sites
avoid_chr = {}
avoid_chr = array(subset(ddply(subset(df.somy.counts, somy == "trisomic" & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)}), V1 < (mean(ddply(subset(df.somy.counts, !(somy %in% c("monosomic", "disomic")) & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)})$V1)*0.1))$chromosome)
avoid_chr = array(subset(ddply(subset(df.somy.counts, somy == "trisomic" & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)}), V1 < (mean(ddply(subset(df.somy.counts, !(somy %in% c("monosomic", "disomic")) & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)})$V1)*0.1))$chromosome)
avoid_chr = unique(c(avoid_chr, array(subset(ddply(df.allele.std, .(chromosome), function(d){dim(d)}), V1 < (mean(ddply(df.allele.std, .(chromosome), function(d){dim(d)})$V1) * 0.3))$chromosome)))
avoid_chr = c(29,avoid_chr)
select_chr1 = c(22, 13, 15, 14, 20, 23, 9, 26)
select_chr2 = c(1, 2, 4, 6, 16, 12, 35, 8, 11, 33, 5)
select_chr3 = c(5,33,3,10,11,8,35,12,18,6,4,16,2, 1)
select_chr = c(1:36)[!(1:36 %in% avoid_chr)]


lm1.sum <- summary(lm(value~surface_ratio,data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness" & (chromosome %in% select_chr1))))
r2.1 <- lm1.sum$r.squared
f1 <- lm1.sum$fstatistic
p1 <- pf(f1[1],f1[2],f1[3],lower.tail=F)
lm2.sum <- summary(lm(value~surface_ratio,data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness" & (chromosome %in% select_chr2))))
r2.2 <- lm2.sum$r.squared
f2 <- lm2.sum$fstatistic
p2 <- pf(f2[1],f2[2],f2[3],lower.tail=F)
lm3.sum <- summary(lm(value~surface_ratio,data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, variable == "std_skewness" & (chromosome %in% select_chr3))))
r3.1 <- lm3.sum$r.squared
f3 <- lm3.sum$fstatistic
p3 <- pf(f3[1],f3[2],f3[3],lower.tail=F)


new_f1c <- ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" ), size = 10, aes(colour = anefreq)) +
geom_text(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" )) +
geom_label(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% select_chr3) & variable == "std_skewness" )) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + 
scale_colour_gradient(low = "white", high = "#EE3B3B", guide= FALSE) + labs(x= "Surface ratio", y = "Skewness", fill = "Trisomy frequency" ) +
theme_bw(base_size = 12, base_family = "ArialMT") + theme(text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

stop("Trolololo")
# This one is selected
svg(filename = "~/tmp/surfaceratio_vs_std_allchr.svg", width = 10, height = 10)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" ), size = 10, aes(colour = anefreq)) +
geom_text(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" )) +
geom_label(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% select_chr3) & variable == "std_skewness" )) +
geom_smooth(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" ), method = "lm", se = TRUE, colour = "blue") +
annotate("text", label = paste("P-value =", format(signif(p3), scientific = TRUE), sep = " "), x = 3, y = 0.7, size = 4, colour = "blue") +
annotate("text", label = paste("R^2 ==", signif(r3.1, digits = 4), sep = " "), x = 3, y = 0.65, size = 4, colour = "blue", parse=TRUE) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + 
scale_colour_gradient(low = "white", high = "#EE3B3B", guide= FALSE) + labs(x= "Surface ratio", y = "Skewness", fill = "Trisomy frequency" ) +
theme_classic(base_size = 12, base_family = "ArialMT")
dev.off()

### Randomization
if(FALSE){

ordered_subset <- subset(df.chr_surfaceratio_stdmeasures_anefreq.m[order(df.chr_surfaceratio_stdmeasures_anefreq.m$chromosome),], variable == "std_skewness")
value_array = array(ordered_subset$value)
surface_ratio_array = array(ordered_subset$surface_ratio)

outputVector <- rep(0,1000000)
x <- c(1:36L)
system.time(
	for(i in 1:1000000) {
		s <- x[.Internal(sample(36L, 14L, FALSE, NULL))]
		outputVector[i] <- summary(lm(value_array[s]~surface_ratio_array[s]))$adj.r.squared
	}
)
sum(outputVector > r3.1) / 1000000

system.time(for(i in 1:1000000) x[.Internal(sample(36L, 14L, FALSE, NULL))])
system.time(lapply(1:1000000, function(i){x[.Internal(sample(36L, 14L, FALSE, NULL))]}))

sum(r3.1 < replicate( 1000, my.fun() )) / 1000
sum(r3.1 < replicate( 10000, my.fun() )) / 10000
sum(r3.1 < replicate( 100000, my.fun() )) / 100000
sum(r3.1 < replicate( 1000000, my.fun() )) / 1000000
}

# Extra figures
if(FALSE) {
postscript("~/tmp/measures_std.eps", width = 20, height = 15, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")
ggplot( df.chr_surfaceratio_stdmeasures_anefreq.m, aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + scale_fill_gradient(low = "white", high = "#EE3B3B") + facet_wrap(~variable, scales="free_y") + labs(x= "Surface ratio", y = "Median standard deviation", fill = "Aneuploidy frquency" )
dev.off()




postscript("~/tmp/ecdf.eps", width = 4, height = 20, horizontal = FALSE, onefile = FALSE, paper = "special", 
	colormodel = "rgb", family = "ArialMT")

                ggplot(subset(df.het.by.somy.plot, somy %in% c("std")), aes(frequency)) + 
                stat_ecdf(geom = "step") +
		scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
		facet_grid(chromosome~somy, scales="free_y") + 
		labs(x = "std") +
		theme_bw(base_size = 8, base_family = "ArialMT") +
		theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())
		
dev.off()


postscript("~/tmp/std_nsamples.eps", width = 12, height = 12, horizontal = FALSE, onefile = FALSE, paper = "special",
    colormodel = "rgb", family = "ArialMT")
    ggplot(df.std.nsamples, aes(samples)) + 
    geom_histogram() +
    facet_wrap(~chromosome, scales="free")
dev.off()


postscript("~/tmp/std_dist_filtered_morethan2samples.eps", width = 15, height = 20, horizontal = FALSE, onefile = FALSE, paper = "special",
    colormodel = "rgb", family = "ArialMT")
    
    ggplot(subset(df.std.nsamples, samples > 2), aes(std)) + 
    geom_histogram() +
    facet_wrap(~chromosome, scales="free")
    
dev.off()

postscript("~/tmp/std_dist_filtered_morethan3samples.eps", width = 15, height = 20, horizontal = FALSE, onefile = FALSE, paper = "special",
    colormodel = "rgb", family = "ArialMT")
    
    ggplot(subset(df.std.nsamples, samples > 3), aes(std)) + 
    geom_histogram() +
    facet_wrap(~chromosome, scales="free")
    
dev.off()

postscript("~/tmp/std_dist_unfiltered_nsamples.eps", width = 15, height = 20, horizontal = FALSE, onefile = FALSE, paper = "special",
    colormodel = "rgb", family = "ArialMT")
    
    ggplot(df.std.nsamples, aes(std)) + 
    geom_histogram() +
    facet_wrap(~chromosome, scales="free")
    
dev.off()
}

# All of thempng("qc-dispersions.png", 1000, 1000, pointsize=20)
png("~/tmp/all_chr_surfacevsskew.png", 1000, 1000, pointsize=20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m,  variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + 
scale_fill_gradient(low = "white", high = "#EE3B3B") + labs(x= "Surface ratio", y = "Skewness", fill = "Aneuploidy frquency" )
dev.off()

png("~/tmp/nofewtrisomics_chr_surfacevsskew.png", 1000, 10000, pointsize=20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + 
scale_fill_gradient(low = "white", high = "#EE3B3B") + labs(x= "Surface ratio", y = "Skewness", fill = "Aneuploidy frquency" )
dev.off()

png("~/tmp/nofewtrisomics_chr_surfacevsskew.png", 1000, 10000, pointsize=20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + 
scale_fill_gradient(low = "white", high = "#EE3B3B") + labs(x= "Surface ratio", y = "Skewness", fill = "Aneuploidy frquency" )
dev.off()


svg(filename = "~/tmp/diffpeaks_vs_std_filters_test.svg", width = 20, height = 20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr)), 
	aes(diff_peaks, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + 
geom_smooth(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness" & (chromosome %in% select_chr1)), method = "lm", se = TRUE,  colour = "blue") +
geom_smooth(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness" & (chromosome %in% select_chr2)), method = "lm", se = TRUE, colour = "blue") +
#annotate("text", label = paste("P-value =", signif(p1, digits = 4), sep = " "), x = 3.2, y = 1.25, size = 8, colour = "black") +
#annotate("text", label = paste("R^2 ==", signif(r2.1, digits = 4), sep = " "), x = 3.2, y = 1.2, size = 8, colour = "black", parse=TRUE) +
#annotate("text", label = paste("P-value =", signif(p2, digits = 4), sep = " "), x = 3, y = 0.7, size = 8, colour = "black") +
#annotate("text", label = paste("R^2 ==", signif(r2.2, digits = 4), sep = " "), x = 3, y = 0.65, size = 8, colour = "black", parse=TRUE) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + labs(x= "Surface ratio", y = "Skewness", fill = "Aneuploidy frquency" ) +
facet_wrap(~variable, scales = "free")
dev.off()

svg(filename = "~/tmp/surfaceratio_vs_std_allchr_nofill.svg", width = 20, height = 20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome)) + 
geom_point() + geom_label() + 
geom_smooth(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, (chromosome %in% select_chr3) & variable == "std_skewness" ), method = "lm", se = TRUE, colour = "blue") +
annotate("text", label = paste("P-value =", format(signif(p3), scientific = TRUE), sep = " "), x = 2.94, y = 0.55, size = 8, colour = "black") +
annotate("text", label = paste("R^2 ==", signif(r3.1, digits = 4), sep = " "), x = 2.94, y = 0.5, size = 8, colour = "black", parse=TRUE)
dev.off()

# Cedric proposed
svg(filename = "~/tmp/surfaceratio_vs_std_filters_test3.svg", width = 20, height = 20)
ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point() + geom_label() + 
geom_smooth(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m, !(chromosome %in% avoid_chr) & variable == "std_skewness"), method = "lm", se = TRUE,  colour = "blue") +
annotate("text", label = paste("P-value =", signif(p3, digits = 4), sep = " "), x = 3.2, y = 1.25, size = 8, colour = "black") +
annotate("text", label = paste("R^2 ==", signif(r3.1, digits = 4), sep = " "), x = 3.2, y = 1.2, size = 8, colour = "black", parse=TRUE) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + labs(x= "Surface ratio", y = "Skewness", fill = "Aneuploidy frquency" )
#facet_wrap(~variable, scales = "free")
dev.off()


rownames(df.chr_surfaceratio_stdmeasures_anefreq) <- df.chr_surfaceratio_stdmeasures_anefreq$chromosome
sub2cluster = df.chr_surfaceratio_stdmeasures_anefreq[,c("surface_ratio", "std_skewness", "anefreq")]

# Typical clustering 
autoplot(kmeans(sub2cluster[, c("surface_ratio", "std_skewness")], 2), data = sub2cluster, label = TRUE, label.size = 3)
autoplot(clara(sub2cluster[, c("surface_ratio", "std_skewness")], 2), data = sub2cluster, label = TRUE, label.size = 3)
autoplot(fanny(sub2cluster[, c("surface_ratio", "std_skewness")], 2), frame = TRUE, data = sub2cluster, label = TRUE, label.size = 3)
autoplot(pam(sub2cluster[, c("surface_ratio", "std_skewness")], 2), data = sub2cluster, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm')

# Local Fisher Discriminant Analysis (LFDA)
model <- lfda(sub2cluster[, c("surface_ratio", "std_skewness")], sub2cluster[, c("surface_ratio", "std_skewness")], r = 2, metric="plain")
autoplot(model, data = sub2cluster, frame = TRUE, frame.colour = 'Species')

# Semi-supervised Local Fisher Discriminant Analysis (SELF)
model <- self(sub2cluster[, c("surface_ratio", "std_skewness")], sub2cluster[, c("surface_ratio", "std_skewness")], beta = 0.1, r = 3, metric="plain")
autoplot(model, data = iris, frame = TRUE, frame.colour = 'Species')


# Make figure S1 allele profiles just divided in two panels given by the selection of chromosomes we are showing here
# select_chr2 <- selected chromosomes for regression 1
# select_chr1 <- selected chromosomes for regression 2

# Sorting ordering for chromosomes
df.het.by.somy.plot.sorted <- df.het.by.somy.plot
ordered_levels = c(array(df.merged$chromosome[order(df.merged$surface_ratio, decreasing=TRUE)]) ,levels(df.het.by.somy.plot$chromosome)[!(levels(df.het.by.somy.plot$chromosome) %in% unique(array(df.merged$chromosome)))])
df.het.by.somy.plot.sorted$chromosome <- factor(df.het.by.somy.plot$chromosome, levels = ordered_levels, ordered = TRUE)

df.het.by.somy.plot.sorted.std <- subset(df.het.by.somy.plot, !(chromosome %in% c(select_chr3, "MC", "MI")) & somy == "std")
df.het.by.somy.plot.sorted.std$chromosome = factor(df.het.by.somy.plot.sorted.std$chromosome, levels = ordered_levels [which(!(ordered_levels %in% c(select_chr3, "MC", "MI")))], ordered = TRUE)


fs1A.1 <- ggplot(subset(df.het.by.somy.plot.sorted, chromosome %in% c(select_chr3) & somy %in% c("disomic", "trisomic")), aes(frequency)) + 
	geom_density(adjust=1.5) +
	facet_grid(chromosome~somy, scales="free") + 
	theme_bw(base_size = 8, base_family = "ArialMT") +
	theme(strip.text.y = element_blank(), strip.background = element_blank())

fs1A.2 <- ggplot(subset(df.het.by.somy.plot.sorted, chromosome %in% c(select_chr3) & somy == "std"), aes(frequency)) + 
	geom_histogram(binwidth=0.01) + 
	scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
	facet_grid(chromosome~., scales="free_y") + 
	labs(x = "std") +
	theme_bw(base_size = 8, base_family = "ArialMT") +
	theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())

fs1B.1 <- ggplot(subset(df.het.by.somy.plot.sorted, !(chromosome %in% c(select_chr3)) & somy %in% c("disomic", "trisomic")), aes(frequency)) + 
	geom_density(adjust=1.5) +
	facet_grid(chromosome~somy, scales="free") + 
	theme_bw(base_size = 8, base_family = "ArialMT") +
	theme(strip.text.y = element_blank(), strip.background = element_blank())

fs1B.2 <- ggplot(subset(df.het.by.somy.plot.sorted.std, !(chromosome %in% c(select_chr3))), aes(frequency)) + 
	geom_histogram(binwidth=0.01) + 
	scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
	facet_grid(chromosome~., scales="free_y", drop = FALSE) + 
	labs(x = "std") +
	theme_bw(base_size = 8, base_family = "ArialMT") +
	theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())

#postscript("~/tmp/newfig1.eps", width = 8, height = 18, horizontal = FALSE, onefile = FALSE, paper = "special", 
#	colormodel = "rgb", family = "ArialMT")

	
	
# To have this figures correctly processed, I need to use the same cutoffs
png("~/tmp/newfigs1.png", width = 600, height = 1800)
	ggdraw() +
  	draw_plot(fs1B.1, 0, 0, (1/3)*2, (2/3)-0.06) +
  	draw_plot(fs1B.2, 0 + (1/3)*2, 0.002, (1/3), (2/3) - 0.07) +
  	draw_plot(fs1A.1, 0, ((1/3)*2)-0.05, (1/3)*2, (1/3) + 0.05) +
  	draw_plot(fs1A.2, 0 + (1/3)*2, (((1/3)*2) + 0.002)-0.05, (1/3), ((1/3)-0.01) + 0.05) +
  	draw_plot_label(c("a","b"), c(0,0), c(1,((1/3)*2)-0.05), size = 15)
dev.off()



# My own comparison between stds

##merged
read.table("merged_before_frequencies/all.piles.depths.txt") -> df.all2merge.piles
colnames(df.all2merge.piles) <- c("chromosome", "start", "allelic_depth", "allele", "name", "depth")

df.all2merge.piles.somy <- merge(df.all2merge.piles, df.het.by.somy[,c("chromosome", "start", "name", "somy")], by = c("chromosome", "start", "name"))

df.merged.frequencies.piles <- ddply(df.all2merge.piles.somy, .(chromosome, start, allele, somy), function(d){ c(sum(d$allelic_depth)/sum(d$depth), dim(d)[1]) })
df.merged.frequencies.piles$chromosome <- factor(array(df.merged.frequencies.piles$chromosome), levels = unique(mixedsort(array(df.merged.frequencies.piles$chromosome))), ordered = TRUE)


df.all2merge.piles.somy.resume <- ddply(df.all2merge.piles.somy, .(chromosome, start, allele), function(d) { 
		print(c(unique(d$chromosome), unique(d$start), unique(d$allele)))
		d.s <- 0
		t.s <- 0
		if(dim(subset(d, somy == "disomic"))[1] >0) {
			d.s <- subset(d, somy == "disomic")$allelic_depth/sum(subset(d, somy == "disomic")$depth)
		}
		if(dim(subset(d, somy == "trisomic"))[1] >0) {
			t.s <- subset(d, somy == "trisomic")$allelic_depth/sum(subset(d, somy == "trisomic")$depth)
		}		

		c(	dim(subset(d, somy == "disomic"))[1],
			dim(subset(d, somy == "trisomic"))[1],
			d.s,
			t.s
		)
	}
)

ddply(subset(df.all2merge.piles.somy.resume, V1 > 1 & V2 > 1), .(chromosome), function(d) { length(unique(d$start)) } )



png("~/tmp/assimilated_profiles0.png", height = 4000, width = 400, pointsize=20)
	ggplot(subset(df.merged.frequencies.piles, somy %in% c("disomic", "trisomic")), aes(V1)) + 
	geom_histogram() + 
	geom_vline(aes(xintercept = 0.25)) + 
	geom_vline(aes(xintercept = 0.33, color = "red")) + 
	geom_vline(aes(xintercept = 0.5)) + 
	geom_vline(aes(xintercept = 0.66, color = "red")) + 
	geom_vline(aes(xintercept = 0.75)) + 
	facet_grid(chromosome~somy, scales = "free")
dev.off()


 
# Sort by aggregated
















#################### Full passage density plot


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


svg(filename="~/plotsave/passage_density_full.svg", height = 20, width = 7)
ggplot(subset(test.p, !Chromosome %in% c("MC", "MI")), aes(x = Frequency)) + geom_density() + facet_grid(Chromosome~Sample) + labs(x = "Frequency", y = "Density", title="") + theme_bw(base_size = 12, base_family = "Helvetica")
dev.off()
