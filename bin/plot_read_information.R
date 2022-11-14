# script to plot read length distribution and other stats about ccs algorithm

# Libraries
library(data.table)
library(ggplot2)
library(stringr)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

# Input files
inp_stats = fread(args[1], h=T, stringsAsFactors=F, sep=",")
#inp_stats = fread("/project/holstegelab/Share/pacbio/data_processed/ad_centenarians/m64367e_221018_215644.ccs.reads_summary.txt", h=T, stringsAsFactors=F, sep=",")
outname = str_replace_all(args[1], '.txt', '.png')

# Main
colnames(inp_stats) = c('Index', 'Read_Name', 'Length', 'Pass', 'Quality', 'Type')
inp_stats$Type[which(inp_stats$Type == 'ccs')] = 'HiFi'
inp_stats$Type[which(inp_stats$Type == 'nonccs')] = 'non-HiFi'
plt1 = ggplot(data = inp_stats, aes(x = Length, fill = Type)) + geom_density(alpha = 0.65) + theme(legend.position = 'top') + xlab('Read length') + ylab('Density') + 
    ggtitle('Read length distribution in HiFi and non-HiFi reads: Density')
plt2 = ggplot(data = inp_stats, aes(x = Type, y = Length, fill = Type)) + geom_violin() + geom_boxplot(width=0.1, color="black", alpha=0.2) + theme(legend.position = 'top') + xlab('Read length') + ylab('Density') + 
    ggtitle('Read length distribution in HiFi and non-HiFi reads: Violin')
plt3 = ggplot(data = inp_stats, aes(x = Pass, y = Quality, color = Type)) + geom_point(stat = 'identity', size = 3, alpha = 0.35) + theme(legend.position = "top") + 
    geom_hline(yintercept = 0.99, linetype="dashed", color = "black") + geom_vline(xintercept = 3, linetype="dashed", color = "black") +
    ggtitle('Number of passes and read quality in Hifi and non-HiFi reads')
combined = ggarrange(plotlist = list(plt1, plt2, plt3), nrow = 3, ncol = 1)
png(outname, height = 18, width = 6, units = 'in', res = 600)
combined
dev.off()