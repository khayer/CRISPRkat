library(reshape)
library(ggplot2)
require(cowplot)
# plot counts 

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

#counts_file = args[1]
#print(counts_file)
#output_file = args[2]
#g1_name = args[3]
#g2_name = args[4]
#sample_names1 = unlist(strsplit(args[5],","))
#sample_names2 = unlist(strsplit(args[6],","))



in_file = args[1] #"~/Dropbox (Work)/files/Abby/GeCKO_Sep_lab_meeting/code/data/median_norm_csgrnas.count_normalized.txt"
title = args[2] #"not normalized"
essential_genes_file = args[3] # "~/Documents/Matt/Abby/GeCKO/essential_genes/core-essential-genes-sym_HGNCID"
prefix = args[4] # "~/Desktop/test"

all_counts = read.csv(in_file, sep ="\t" , stringsAsFactors = F)
head(all_counts)
essential_genes = read.csv(essential_genes_file, sep ="\t" , header = F)$V1
head(essential_genes)
#cummulative distribution
tmp_df <- melt(all_counts[,-(1:2)]+1)
head(tmp_df)
cummalative <- ggplot(tmp_df, aes(value, colour = variable)) + 
  stat_ecdf() + 
  scale_x_continuous(trans = 'log2') +
  ggtitle(title) +
  theme_bw() +
  xlab("Log2 counts") +
  ylab("Coverage of sgRNA library\n(cumulative frequency)")
cummalative

ggsave(paste(prefix, "cumalitve.pdf",sep = "_"), cummalative, width = 6, height = 5 )

tmp_df = melt(all_counts)
tmp_df$kind = "other"
if (dim(tmp_df[tmp_df$Gene %in% essential_genes,])[1] > 0 ) {
  tmp_df[tmp_df$Gene %in% essential_genes,]$kind = "essential"
}
print("BLABLA")
print(dim(tmp_df[grep("CTRL0", tmp_df$sgRNA),]))
if (dim(tmp_df[grep("NTCONTROL", tmp_df$Gene),])[1] > 0) {
  tmp_df[grep("NTCONTROL", tmp_df$Gene),]$kind = "non targeting"
} else if (dim(tmp_df[grep("NonTar", tmp_df$Gene),])[1] > 0) {
  tmp_df[grep("NonTar", tmp_df$Gene),]$kind = "non targeting"
} else if (dim(tmp_df[grep("CTRL0", tmp_df$sgRNA),])[1] > 0) {
  tmp_df[grep("CTRL0", tmp_df$sgRNA),]$kind = "non targeting"
}

box_plot = ggplot(tmp_df, aes(kind,log2(value+1), fill = variable)) +
  geom_boxplot() + theme_bw() +
  ggtitle(title) +
  ylab("sgRNA representation\nlog2(count+1)") +
  xlab("kind of sgRNA")
box_plot

ggsave(paste(prefix, "boxplot.pdf",sep = "_"), box_plot, width = 6, height = 5 )

plot_grid(cummalative, box_plot, labels = "AUTO",align = 'h', label_size = 12)
ggsave(paste(prefix, "combined.pdf",sep = "_"), width = 10, height = 5 )
