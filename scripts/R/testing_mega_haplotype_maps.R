

plots <- list()
i = 1
	#			 1    3    5    7    9     11    13    15    17    19    21    23			
	for (c in c("3", "5", "8", "9", "14", "15", "16", "20", "23", "26", "31", "32")) {  
		subsetmelt = subset(df.sites.melt, chr == c)
		ntimes = length(levels(subset(df.sites.melt, chr == c)$variable))
		nmaxpos = (dim(subset(df.sites.melt, chr == c))[1]/ntimes)
		nbreaks = c(1,round(nmaxpos/2),nmaxpos)
		subsetmelt$npos = rep(c(1:nmaxpos), times = ntimes)
		plots[[i]] <- ggplot(subsetmelt, aes(y=variable, x=npos)) + 
		geom_tile(aes(fill=factor(value), height = .7)) + 
		scale_fill_manual(values=pcolors, labels=c("A","C","G","T"), breaks = c(0:3)) +
		scale_x_continuous(expand = c(0, 0), breaks=nbreaks) + 
		facet_grid(.~chr, scales = "free") +
		#labs(title = "Genome haplotype map", y = "Clone", x = "Genome position", fill = "Allele") +
		labs(title = "", y = "", x = "", fill = "") +
		theme_bw(base_size = 16, base_family = "ArialMT") +
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

		plots[[i+1]] <- ggplot(subset(df.pclones.tiles.count.somy, chromosome == c), aes(sample, count)) + 
		geom_bar(stat = "identity", width=.7, fill = 'white', colour = 'black') + 
		#scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		coord_flip() + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="", fill = "somy") + 
		theme_bw(base_size = 16, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 16, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(), legend.position = "none",
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank())
        	
                plots[[i+2]] <- ggplot(subset(df.pclones.piles.medians, chromosome == c), aes(frequency, fill = somy)) + 
		geom_density() + scale_fill_gradient( limits=c(1.6, 4.1) , breaks = c(2, 3, 4), labels = c("2", "3", "4"), low = "red", high = "green") + 
		scale_x_continuous(breaks = c(0.25,0.5,0.75)) + 
		facet_wrap(~sample, ncol = 1 ) + 
		#labs(x = "", y="", fill = "Somy") + 
		labs(x = "", y="", fill = "somy") + 
		theme_bw(base_size = 16, base_family = "ArialMT") +
		theme(plot.background = element_rect(fill = "black") , text=element_text( size = 16, family="ArialMT"), 
			plot.margin = unit(c(0,0,-1,-1), "cm"), panel.margin = unit(0, "lines"), 
			axis.text.x = element_text(angle = 90), strip.background = element_blank(), 
			strip.text.x = element_blank(),
        	axis.text.y=element_blank(),
        	axis.ticks.y=element_blank(), legend.position = "none")
        	
		i = i + 3
	}

	### TEmporal solo para plotear ms rapido
	png("~/tmp/test_new_megahaplomap_submit.png", 1000, 1000, pointsize=20)
	ggdraw() +

	draw_plot(plots[[4]],  0.02 + wc3, 	yp11, wc1, hc1+0.025) + 
	draw_plot(plots[[10]], 0.335  + wc3, yp11, wc1, hc1+0.025) + 
	draw_plot(plots[[13]], 0.665+ wc3, yp11, wc1, hc1+0.025) + 

	draw_plot(plots[[6]], 0.03, yp21, wc3, hc2) +	
	draw_plot(plots[[12]], 0.34 , yp21, wc3, hc2) +
	draw_plot(plots[[15]], 0.67 , yp21, wc3, hc2) +


	draw_plot(plots[[5]], wc3 + 0.2 + 0.025, yp21, wc2, hc2) +	
	draw_plot(plots[[11]], wc3 + 0.2 + 0.335, yp21, wc2, hc2) +
	draw_plot(plots[[14]], wc3 + 0.2 + 0.665, yp21, wc2, hc2) +

##
	draw_plot(plots[[16]],  0.025 + wc3, 	yp12, wc1, hc1+0.015) + 
	draw_plot(plots[[22]], 0.335  + wc3, 	yp12, wc1, hc1+0.015) + 
	draw_plot(plots[[25]], 0.665+ wc3, 		yp12, wc1, hc1+0.015) + 

	draw_plot(plots[[18]], 0.03, yp22, wc3, hc2) +	
	draw_plot(plots[[24]], 0.34 , yp22, wc3, hc2) +
	draw_plot(plots[[27]], 0.67 , yp22, wc3, hc2) +


	draw_plot(plots[[17]], wc3 + 0.2 + 0.025, yp22, wc2, hc2) +	
	draw_plot(plots[[23]], wc3 + 0.2 + 0.335, yp22, wc2, hc2) +
	draw_plot(plots[[26]], wc3 + 0.2 + 0.665, yp22, wc2, hc2) +

##
	draw_plot(plots[[28]],  0.025 + wc3, 	yp13, wc1, hc1+0.015) + 
	draw_plot(plots[[31]], 0.335  + wc3, 	yp13, wc1, hc1+0.015) + 

	draw_plot(plots[[30]], 0.03, yp23, wc3, hc2) +	
	draw_plot(plots[[33]], 0.34 , yp23, wc3, hc2) +


	draw_plot(plots[[29]], wc3 + 0.2 + 0.025, yp23, wc2, hc2) +	
	draw_plot(plots[[32]], wc3 + 0.2 + 0.335, yp23, wc2, hc2) +

	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +

	draw_plot(f4b, 0, 0, 1, 0.26) +
        draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.26), size = 16)
        
dev.off()

#### TMP final ###

	
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

  leg1 <- g_legend(tmplot1)
  leg2 <- g_legend(tmplot2)


xp1 = 0.25
yp1 = 0.265
yp2 = 0.272
wc1 = 0.18
diffp = 0.03
wc2 = 0.035
wc3 = 0.035
hc1 = 0.2
hc2 = 0.17
hc3 = 0.15

twc = wc1+wc2+wc3


# Which chromsomes want to show
#  1    3    5    7     9    11    13    15    17    19    21    23   
#  1    4    7    10   13    16    19    22    25    28    31    34 
# "3", "5", "8", "9", "14", "15", "16", "20", "23", "26", "31", "32"
#

yp11 = (yp1*3) - 0.04
yp21 = (yp2*3) - 0.046

yp12 = (yp1*2) - 0.01
yp22 = (yp2*2) - 0.01

yp13 = (yp1) + 0.04
yp23 = (yp2) + 0.045



png("~/tmp/test_new_megahaplomap_b_submit.png", 1000, 200, pointsize=20)
f4b <- ggplot(subset(test.p, Chromosome %in% c("5", "9", "14", "15", "20", "23", "26", "31")), aes(x = Frequency, colour = Sample)) + 
	geom_density() + 
	facet_wrap(~Chromosome, ncol=8) + 
	labs(x = "Frequency", y = "Density", title="") + 
	theme_bw(base_size = 12, base_family = "ArialMT") + 
	science_theme2 +
	theme(legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') +
	guides(colour=guide_legend(nrow=1,byrow=FALSE), size = guide_legend(order = 3))
#	theme(
#		legend.margin = unit(c(-0.6,0,0,0), "lines"),
#    	#axis.text.x=element_blank(), axis.ticks.x=element_blank(),  
#	    legend.title  = element_text(face = 'bold'), legend.direction = 'horizontal', legend.position = 'top') + 
#		guides(fill=guide_legend(nrow=1,byrow=FALSE), linetype = guide_legend(nc = 1,keywidth=20))
dev.off()

png("~/tmp/test_new_megahaplomap_submit.png", 1000, 1000, pointsize=20)
	ggdraw() +
	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +

	draw_plot(plots[[4]],  0.02 + wc3, 	yp11, wc1, hc1+0.015) + 
	draw_plot(plots[[10]], 0.335  + wc3, yp11, wc1, hc1+0.015) + 
	draw_plot(plots[[13]], 0.665+ wc3, yp11, wc1, hc1+0.015) + 

	draw_plot(plots[[6]], 0.03, yp21, wc3, hc2) +	
	draw_plot(plots[[12]], 0.34 , yp21, wc3, hc2) +
	draw_plot(plots[[15]], 0.67 , yp21, wc3, hc2) +


	draw_plot(plots[[5]], wc3 + 0.2 + 0.025, yp21, wc2, hc2) +	
	draw_plot(plots[[11]], wc3 + 0.2 + 0.335, yp21, wc2, hc2) +
	draw_plot(plots[[14]], wc3 + 0.2 + 0.665, yp21, wc2, hc2) +

##
	draw_plot(plots[[16]],  0.025 + wc3, 	yp12, wc1, hc1+0.015) + 
	draw_plot(plots[[22]], 0.335  + wc3, 	yp12, wc1, hc1+0.015) + 
	draw_plot(plots[[25]], 0.665+ wc3, 		yp12, wc1, hc1+0.015) + 

	draw_plot(plots[[18]], 0.03, yp22, wc3, hc2) +	
	draw_plot(plots[[24]], 0.34 , yp22, wc3, hc2) +
	draw_plot(plots[[27]], 0.67 , yp22, wc3, hc2) +


	draw_plot(plots[[17]], wc3 + 0.2 + 0.025, yp22, wc2, hc2) +	
	draw_plot(plots[[23]], wc3 + 0.2 + 0.335, yp22, wc2, hc2) +
	draw_plot(plots[[26]], wc3 + 0.2 + 0.665, yp22, wc2, hc2) +

##
	draw_plot(plots[[28]],  0.025 + wc3, 	yp13, wc1, hc1+0.015) + 
	draw_plot(plots[[31]], 0.335  + wc3, 	yp13, wc1, hc1+0.015) + 

	draw_plot(plots[[30]], 0.03, yp23, wc3, hc2) +	
	draw_plot(plots[[33]], 0.34 , yp23, wc3, hc2) +


	draw_plot(plots[[29]], wc3 + 0.2 + 0.025, yp23, wc2, hc2) +	
	draw_plot(plots[[32]], wc3 + 0.2 + 0.335, yp23, wc2, hc2) +


	draw_plot(f4b, 0, 0, 1, 0.26) +
        draw_plot_label(c("a", "b"), c(0, 0), c(1, 0.26), size = 16)
        
dev.off()


yp14 = (0) + 0.04
yp24 = (0) + 0.052
# Which chromsomes want to show
#  1    3    5    7     9    11    13    15    17    19    21    23   
#  1    4    7    10   13    16    19    22    25    28    31    34 
# "3", "5", "8", "9", "14", "15", "16", "20", "23", "26", "31", "32"
png("~/tmp/test_new_megahaplomap_all.png", 1000, 1000, pointsize=20)
	ggdraw() +
	draw_grob(leg1, 0.25, 0.95, 0.25, 0.05) + 
	draw_grob(leg2, 0.5, 0.95, 0.25, 0.05) +

	draw_plot(plots[[1]],  0.02 + wc3, 	yp11, wc1, hc1) + 
	draw_plot(plots[[4]], 0.335  + wc3, yp11, wc1, hc1) + 
	draw_plot(plots[[7]], 0.665+ wc3, yp11, wc1, hc1) +

	draw_plot(plots[[3]], 0.03, yp21, wc3, hc2) +	
	draw_plot(plots[[6]], 0.34 , yp21, wc3, hc2) +
	draw_plot(plots[[9]], 0.67 , yp21, wc3, hc2) +


	draw_plot(plots[[2]], wc3 + 0.2 + 0.025, yp21, wc2, hc2) +	
	draw_plot(plots[[5]], wc3 + 0.2 + 0.335, yp21, wc2, hc2) +
	draw_plot(plots[[8]], wc3 + 0.2 + 0.665, yp21, wc2, hc2) +

##
	draw_plot(plots[[10]],  0.025 + wc3, 	yp12, wc1, hc1) + 
	draw_plot(plots[[13]], 0.335  + wc3, 	yp12, wc1, hc1) + 
	draw_plot(plots[[16]], 0.665+ wc3, 		yp12, wc1, hc1) +

	draw_plot(plots[[12]], 0.03, yp22, wc3, hc2) +	
	draw_plot(plots[[15]], 0.34 , yp22, wc3, hc2) +
	draw_plot(plots[[18]], 0.67 , yp22, wc3, hc2) +


	draw_plot(plots[[11]], wc3 + 0.2 + 0.025, yp22, wc2, hc2) +	
	draw_plot(plots[[14]], wc3 + 0.2 + 0.335, yp22, wc2, hc2) +
	draw_plot(plots[[17]], wc3 + 0.2 + 0.665, yp22, wc2, hc2) +

##
	draw_plot(plots[[19]],  0.025 + wc3, 	yp13, wc1, hc1) + 
	draw_plot(plots[[22]], 0.335  + wc3, 	yp13, wc1, hc1) + 
	draw_plot(plots[[25]], 0.665  + wc3, 	yp13, wc1, hc1) + 

	draw_plot(plots[[21]], 0.03, yp23, wc3, hc2) +	
	draw_plot(plots[[24]], 0.34 , yp23, wc3, hc2) +
	draw_plot(plots[[27]], 0.67 , yp23, wc3, hc2) +


	draw_plot(plots[[20]], wc3 + 0.2 + 0.025, yp23, wc2, hc2) +	
	draw_plot(plots[[23]], wc3 + 0.2 + 0.335, yp23, wc2, hc2) +
	draw_plot(plots[[26]], wc3 + 0.2 + 0.665, yp23, wc2, hc2) +

##
	draw_plot(plots[[28]],  0.025 + wc3, 	yp14, wc1, hc1) + 
	draw_plot(plots[[31]], 0.335  + wc3, 	yp14, wc1, hc1) + 
	draw_plot(plots[[34]], 0.665  + wc3, 	yp14, wc1, hc1) + 

	draw_plot(plots[[30]], 0.03, yp24, wc3, hc2) +	
	draw_plot(plots[[33]], 0.34 , yp24, wc3, hc2) +
	draw_plot(plots[[36]], 0.67 , yp24, wc3, hc2) +

	draw_plot(plots[[29]], wc3 + 0.2 + 0.025, yp24, wc2, hc2) +	
	draw_plot(plots[[32]], wc3 + 0.2 + 0.335, yp24, wc2, hc2) +
	draw_plot(plots[[35]], wc3 + 0.2 + 0.665, yp24, wc2, hc2)

dev.off()