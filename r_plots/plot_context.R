library(ggplot2)


# Rscript plotContext.R input_file output_dir
args = commandArgs(trailingOnly=TRUE)
inputFile = args[1]
outputDir = args[2]

# counts inputFile =  "/home/enric_serra/Projects/reportCancer/report_cancer/TEST/snv_context.txt"
mytable = read.table(inputFile, sep="\t", header=FALSE)
humanGenomeFrequency = c(119012680, 67747346, 15146853, 94892673, 107932070, 76434609, 16303729, 103793269, 84736529, 69175754, 14111609,
                          81848639, 115628225, 91105942, 13244230, 130932790, 120766278, 78713545, 107729359, 146205362, 75676004, 99106720,
                          118939635, 117744923, 66253092, 55279316, 88658468, 86608671, 120765137, 117811750, 111982043, 226379708)
counts = as.numeric(mytable[1, 1:96]) / humanGenomeFrequency
counts = counts/sum(counts)*100
countsPositionsStart = 1:96 - 0.9
countsPositionsEnd = 1:96 - 0.1

# Position and tags for the x axis
cytosineTags = rep(c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT"), 3)
timineTags = rep(c("ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT"), 3)
countTags = c(cytosineTags, timineTags)
countTagsPosition = 1:96 - 0.5

# Position and tags for the y axis
yAxis = 0:4 * 5
yAxisText=as.character(yAxis)
yAxis = yAxis * 5

# Get sets of colors, redefine them here, and it will be propagated
plotColors = c("blue", "green", "purple", "orange", "red", "yellow")
gRectColors = c()
for(aColor in plotColors){
  gRectColors = c(gRectColors, rep(aColor, 16))     
}
# Rectangles at the top of the plot
topLevelRectanglesPositionStart = 0:5 * 16 + 0.1
topLevelRectanglesPositionEnd = 1:6 * 16 - 0.1
topLevelTags = c("C > A", "C > G", "C > T" , "T > A", "T > C", "T > G") 
topLevelTagsPosition = 0:5 * 16 + 8 



# Plot
png(paste(outputDir, "/snv_context.png", sep=""), width=1628, height = 1028)

ggplot()+
  
  geom_rect(mapping=aes(xmin=countsPositionsStart, xmax=countsPositionsEnd, ymin=0, ymax=counts), color="white", fill=gRectColors) + 

  geom_rect(mapping=aes(xmin=topLevelRectanglesPositionStart, xmax=topLevelRectanglesPositionEnd, ymin=105, ymax=107), color="white", fill=plotColors) +
  
  geom_text(mapping=aes(x=topLevelTagsPosition, y=112), label=topLevelTags) +
  
  geom_text(mapping=aes(x=countTagsPosition, y=-4), label=countTags, angle=90, size=3) +

  geom_segment(mapping=aes(x=-0.5, xend=-0.5, y=0, yend=100)) +
  
  geom_segment(mapping=aes(x=-1, xend=-0.5, y=yAxis, yend=yAxis)) +
  geom_text(mapping=aes(x=-2, y=yAxis), label=yAxisText) +
  
  xlim(-4, 100) + ylim(-5,120)+
  
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),       
      axis.ticks.y=element_blank(),
      axis.title.x=element_blank(),
       axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      panel.background=element_blank(), 
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank()
  ) +
  geom_text(mapping=aes(x=-4, y=50),label = "Normalised substitutions", angle=90)
  

dev.off()


counts = as.numeric(mytable[1, 1:96])

countsAggregated = c(sum(counts[1:16]))

countsAggregated = c(sum(counts[1:16]), sum(counts[17:32]), sum(counts[33:48]), sum(counts[49:64]), sum(counts[65:80]), sum(counts[81:96]))

xPositionStart = 1:6 -0.9
xPositionEnd = 1:6 -0.1
xTags = c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C")
xTagsPosition = 1:6 - 0.5

nearestThousand = function(number){
  return((number%/%1000)+1)     
}
xTicksPosition = (0:nearestThousand(max(countsAggregated))*1000)
xTicksText= as.character(xTicksPosition)

png(paste(outputDir, "/aggregated_substitutions.png", sep=""), width=1200,height =560)

  ggplot() + geom_rect(mapping=aes(xmin=xPositionStart, xmax=xPositionEnd, ymin=0, ymax=countsAggregated), color ="white",fill= plotColors) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),       
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          panel.background=element_blank(), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank()
    ) +
    geom_segment(mapping=aes(x=-0.05, xend=-0.05, y=0, yend=max(xTicksPosition))) +
    geom_segment(mapping=aes(x=-0.05, xend=-0.15, y=xTicksPosition, yend=xTicksPosition)) +
    geom_text(mapping=aes(x=-0.3,  y=xTicksPosition), label=xTicksText) +
    
    geom_text(mapping=aes(x=xTagsPosition, y=-(max(countsAggregated)* 0.07)), label=xTags, angle=90) +
    xlim(-2.2, 7)
dev.off()