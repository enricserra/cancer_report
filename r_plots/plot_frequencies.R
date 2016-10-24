library(ggplot2)
args = commandArgs(trailingOnly=TRUE)

mytable = read.table(args[1], sep="\t", header=FALSE)

file_name = paste(args[3], "snv_frequencies.png", sep="/")

png(file_name, width=750, height=480)
    qplot(mytable[,1],binwidth = 0.01) + geom_histogram(color="white",fill="red",binwidth=0.01) + theme_classic() +
    scale_x_continuous()+xlab("Allele frequency") + ylab("Counts") + ggtitle("SNV Alelle Frequencies") + xlim(0, 1)
dev.off()

mytable = read.table(args[2], sep="\t", header=FALSE)

file_name = paste(args[3], "indel_frequencies.png", sep="/")

png(file_name,width=750, height=480)
    qplot(mytable[,1],binwidth = 0.01) + geom_histogram(color="white",fill="green",binwidth=0.01) + theme_classic() +
    scale_x_continuous() + xlab("Allele frequency") + ylab("Counts") + ggtitle("Indel Allele Frequencies") + xlim(0,1)
dev.off()
