library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

chromosome_accumulated = c(0, 248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 
                           1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562,
                           2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417, 3088269832 )

isPlottable = function(allele) {
  return((allele=="C" || allele=="T"))
}

reverse = function(allele) {
  if(allele == "G") { return("C") }  
  if(allele=="A") { return("T") }
  if(allele=="T") { return("A") }
  if(allele=="C"){ return("G") }
  return(allele)
}

assignMutation  = function(object) {
  ref = object[1]
  alt = object[2]
  
  if(!isPlottable(ref)) { 
    ref = reverse(ref)
    alt = reverse(alt)
  }
  return(paste(ref, alt, sep=" > "))
}


kataegis_input = read.table(args[1], sep="\t", header=FALSE)

chromosomeDataFrame = as.data.frame(cbind(chromosome_accumulated[seq(from = 2, to = length(chromosome_accumulated), by = 2)], chromosome_accumulated[seq(from = 3, to = length(chromosome_accumulated), by = 2)]))

meanSizeChromosomeDataFrame = as.data.frame(cbind(apply(cbind(chromosome_accumulated[seq(from = 1, to = length(chromosome_accumulated)-1, by =1)], 
                                                              chromosome_accumulated[seq(from = 2, to = length(chromosome_accumulated), by =1)]),1 ,mean),
                                                        paste("chr",c(1:22,"X","Y"),sep="")))

distanceDataFrame = data.frame(kataegis_input)

hotspot_table = read.table(args[2], sep="\t", header=FALSE)

hotspot_dataframe = as.data.frame(hotspot_table)

kataegis_plot = ggplot() + 
  
  geom_rect(data=chromosomeDataFrame, aes(xmax=chromosomeDataFrame[,1], xmin=chromosomeDataFrame[,2],
                                          ymax = rep(max(kataegis_input[,3])*1.15, length(chromosomeDataFrame[,1])),
                                         ymin = rep(0.7,length(chromosomeDataFrame[,2]))), colour="white",fill="#ECF3F6") +
  
  geom_point(data=distanceDataFrame, aes(distanceDataFrame[,2], (distanceDataFrame[,3]),
                                         colour=apply(cbind(as.character(distanceDataFrame[,4]),                
                                                            as.character(distanceDataFrame[,5])),1 , assignMutation))) +

  geom_text(data=meanSizeChromosomeDataFrame,
            mapping=aes(x=as.numeric(as.character(meanSizeChromosomeDataFrame[,1])),y=0.5),angle=90, ,size=3,
            label=meanSizeChromosomeDataFrame[,2], colour="black") +
  
  geom_segment(data=hotspot_dataframe,mapping=aes(x =hotspot_dataframe[,1] , y = 1.5, xend = hotspot_dataframe[,1], yend = 3)) + 
  
  geom_segment(data=hotspot_dataframe,mapping=aes(x =hotspot_dataframe[,1]+10000000 , y = 2.7, xend = hotspot_dataframe[,1], yend = 3)) + 
  
  geom_segment(data=hotspot_dataframe,mapping=aes(x=hotspot_dataframe[,1]-10000000 , y = 2.7, xend = hotspot_dataframe[,1], yend = 3)) + 
  
  geom_segment(mapping=aes(x=0, xend=3088269832, y=0.7, yend=0.7)) +
  
  geom_segment(mapping=aes(x=3088269832, xend=3088269832, y=0.7, yend=0.52)) +
  
  geom_segment(mapping=aes(x=0, xend=0, y=0.7, yend=0.52)) +
  
  geom_segment(mapping=aes(x=-40000000, xend=-40000000, y=1, yend=max(kataegis_input[,3]) * 1.15)) +
  
  geom_segment(mapping=aes(x=-40000000, xend=-65000000, y=1, yend=1)) +
  
  geom_segment(mapping=aes(x=-40000000, xend=-65000000, y=max(kataegis_input[,3])*1.15, yend=max(kataegis_input[,3]) * 1.15))+
  
  xlim(-65000000, chromosome_accumulated[26] * 1.1) + ylim(0.3, max(kataegis_input[,3]) * 1.1) + xlab("") + ylab("") + 
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),

        panel.background=element_blank(), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()
        ) +
  
  scale_y_log10() +
  
  scale_color_manual(values=c("#827BE9","black", "#E47B3E", "pink", "#F0E877", "#A5DF94"),name="Nucleotide\nChange")
  
output_file_name = paste(args[3], "kataegis_with_hotspots.png", sep="/")
ggsave(output_file_name, plot = kataegis_plot, device = "png", path = NULL, scale = 1, width = 20, height = 10, dpi = 300)