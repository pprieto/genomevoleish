	df.allele.std.0 <- read.table("../contrast_alleles/all/all_real_nsamples_results.txt")
	colnames(df.allele.std.0) <- c("chromosome", "start", "allele", "frequency", "samples")



#  Using the allele frequency surface from all sites without unfilter
df.std.maxsamples = ddply(df.allele.std.0, .(chromosome), function(d) { max(d$samples) } )
colnames(df.std.maxsamples) <- c( "chromosome", "max" )
rownames(df.std.maxsamples)  <- array(df.std.maxsamples$chromosome)

df.allele.std.1 <- df.allele.std.0[0,]
# Filter out the others
for(t in c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) {
#for(t in c(0.9)) {
    print(c("Working this shit for threshold ", t))
    df.allele.std.1 <- df.allele.std.0[0,]
    for(c in levels(df.allele.std.0$chromosome)) {
        print(c (c, df.std.maxsamples[c,"max"]*t, dim(subset(df.allele.std.0, chromosome == c & samples > df.std.maxsamples[c,"max"] * t))[1] ))
        df.allele.std.1 <- rbind(df.allele.std.1, subset(df.allele.std.0, chromosome == c & samples > df.std.maxsamples[c,"max"]*t ))
    }

	df.allele.std$somy <- "std"
	df.allele.std.1$somy <- "std"
	df.het.by.somy.plot <- rbind(df.het.by.somy[,c("chromosome", "start", "frequency", "somy")], df.allele.std[,c("chromosome", "start", "frequency", "somy")])
	df.het.by.somy.plot.1 <- rbind(df.het.by.somy[,c("chromosome", "start", "frequency", "somy")], df.allele.std.1[,c("chromosome", "start", "frequency", "somy")])

df.chr_trisomic = ddply(subset(df.somy.counts, (somy == "trisomic") & !is.na(somy)), .(chromosome), function(d) { sum(d$samples)})
colnames(df.chr_trisomic) <- c("chromosome", "trisomic")

df.chr_sites4std = ddply(df.allele.std.1, .(chromosome), function(d){dim(d)[1]})
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

df.meanstd = ddply(subset(df.het.by.somy.plot.1, somy == "std" ), .(chromosome),  function(d) {
	mean(d$frequency, na.rm = TRUE)
	}
)
colnames(df.meanstd) <- c("chromosome" , "mean_std")

df.medstd = ddply(subset(df.het.by.somy.plot.1, somy == "std" ), .(chromosome),  function(d) {
	median(d$frequency, na.rm = TRUE)
	}
)
colnames(df.medstd) <- c("chromosome" , "median_std")

df.stdskewness = ddply(subset(df.het.by.somy.plot.1, somy == "std" ), .(chromosome),  function(d) {
	skewness(d$frequency, na.rm = TRUE)
	}
)
colnames(df.stdskewness) <- c("chromosome" , "std_skewness")

df.stdkurtosis = ddply(subset(df.het.by.somy.plot.1, somy == "std" ), .(chromosome),  function(d) {
	kurtosis(d$frequency, na.rm = TRUE)
	}
)
colnames(df.stdkurtosis) <- c("chromosome" , "std_kurtosis")

colnames(df.aneuploid.field) <- c("chromosome" , "anefreq" )

df.merged <- Reduce(function(x,y) {merge(x,y, by = "chromosome")}, list(df.chr_trisomic, df.chr_sites4std,df.chr_surfaceratio, df.chr_diffpeaks, df.meanstd, df.medstd, df.stdskewness, df.stdkurtosis, df.aneuploid.field))

df.chr_surfaceratio_stdmeasures_anefreq.m = melt(df.merged, id.vars = c("chromosome", "diff_peaks", "surface_ratio", "anefreq"), measure.vars= c("mean_std", "median_std", "std_kurtosis", "std_skewness"))


new_f1d <- ggplot( subset(df.chr_surfaceratio_stdmeasures_anefreq.m, variable == "std_skewness"), 
	aes(surface_ratio, value, label = chromosome, fill = anefreq)) + 
geom_point(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m,  variable == "std_skewness" ), size = 5, aes(colour = anefreq)) +
#geom_text(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m,  variable == "std_skewness" )) +
geom_label_repel(data=subset(df.chr_surfaceratio_stdmeasures_anefreq.m,  variable == "std_skewness" ),
	fontface = 'bold',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50'
	) +
scale_fill_gradient(low = "white", high = "#EE3B3B") + 
scale_colour_gradient(low = "white", high = "#EE3B3B", guide= FALSE) + labs(x= "Surface ratio", y = "Skewness", fill = "Trisomy frequency" ) +
theme_bw(base_size = 12, base_family = "ArialMT") + theme(text = element_text(face = "bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


png(paste("~/tmp/density/trisomics/",t,"/newfigCorr.png", sep = ""), width = 800, height = 800)
	p <- ggdraw() + draw_plot(new_f1c, 0, 0, 1, 1)
	print(p)
	dev.off()

df.het.by.somy.plot.sorted <- df.het.by.somy.plot
ordered_levels = c(array(df.merged$chromosome[order(df.merged$surface_ratio, decreasing=TRUE)]) ,levels(df.het.by.somy.plot$chromosome)[!(levels(df.het.by.somy.plot$chromosome) %in% unique(array(df.merged$chromosome)))])
df.het.by.somy.plot.sorted$chromosome <- factor(df.het.by.somy.plot$chromosome, levels = ordered_levels, ordered = TRUE)

df.het.by.somy.plot.sorted.1 <- df.het.by.somy.plot.1
ordered_levels = c(array(df.merged$chromosome[order(df.merged$surface_ratio, decreasing=TRUE)]) ,levels(df.het.by.somy.plot.1$chromosome)[!(levels(df.het.by.somy.plot.1$chromosome) %in% unique(array(df.merged$chromosome)))])
df.het.by.somy.plot.sorted.1$chromosome <- factor(df.het.by.somy.plot.1$chromosome, levels = ordered_levels, ordered = TRUE)

df.het.by.somy.plot.sorted.1 <- subset(df.het.by.somy.plot.sorted.1, !(chromosome %in% c("MC", "MI")) & somy == "std")
df.het.by.somy.plot.sorted.1$chromosome = factor(df.het.by.somy.plot.sorted.1$chromosome, levels = ordered_levels [which(!(ordered_levels %in% c("MC", "MI")))], ordered = TRUE)

fs1A.1 <- ggplot(subset(df.het.by.somy.plot.sorted, !(chromosome %in% c("MC", "MI")) & somy %in% c("disomic", "trisomic")), aes(frequency)) + 
	geom_density(adjust=1.5) +
	facet_grid(chromosome~somy, scales="free") + 
	theme_bw(base_size = 8, base_family = "ArialMT") #+
	#theme(strip.text.y = element_blank(), strip.background = element_blank())

fs1A.2 <- ggplot(subset(df.het.by.somy.plot.sorted.1, somy == "std"), aes(frequency)) + 
	#geom_histogram(binwidth=0.01) + 
	geom_density() +
	scale_x_continuous(limits = c(0,0.2), expand=c(0,0), breaks = c(0,0.1,0.2)) +
	facet_grid(chromosome~., scales="free_y", drop = FALSE) + 
	labs(x = "std") +
	theme_bw(base_size = 8, base_family = "ArialMT") +
	theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0), strip.text.x = element_blank())

# To have this figures correctly processed, I need to use the same cutoffs
png(paste("~/tmp/density/trisomics/",t,"/newfigs1.png", sep = ""), width = 600, height = 1800)
	p <- ggdraw() +
  	draw_plot(fs1A.1, 0, 0, (1/3)*2, 1) +
  	draw_plot(fs1A.2, 0 + (1/3)*2, 0, (1/3), 1-0.01) + draw_plot_label(c("a","b"), c(0,0), c(0,0), size = 15)
  	print(p)
dev.off()

}


## Filtering also the allele frquency surface
