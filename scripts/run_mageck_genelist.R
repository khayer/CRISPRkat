library(ggplot2)
library(ggrepel)

# Rscript {workflow.basedir}/scripts/run_pbnpa_genelist.R {input[0]} {input[1]} {output[0]} {params.norm_method} {output[1]} {output[2]} 

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

number_of_genes_displayed = 15

control_comp_file = args[1]
treatment_comp_file = args[2]
output_file = args[3]

norm_method = args[4]
negative_selection_plot = args[5]
positive_selection_plot = args[6]
essential = read.csv(args[7],sep="\t", header = F)
all_genes = args[8]

if (all_genes == "True") {
  all_genes = T
} else {
  F
}
gene_set = unlist(strsplit(args[9],","))

# for testing only!!
#control_comp_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/mageck_total_day0_unt_vs_day15_unt.gene_summary.txt"
#treatment_comp_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/mageck_total_day0_unt_vs_day15_low_dox.gene_summary.txt"
#output_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/mageck_test_combined.csv"
#norm_method = "total"
#negative_selection_plot = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/mageck_negative.pdf"
#positive_selection_plot ="/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/mageck_positve.pdf"

ylab_name = gsub(".gene_summary.txt",x = tail(strsplit(treatment_comp_file,"/")[[1]], n=1),replacement = "")
ylab_name = gsub("^mageck_[A-Za-z]*_", ylab_name, replacement = "")

xlab_name = gsub(".gene_summary.txt",x = tail(strsplit(control_comp_file,"/")[[1]], n=1),replacement = "")
xlab_name = gsub("^mageck_[A-Za-z]*_", xlab_name, replacement = "")

control_comp = read.csv(control_comp_file,  stringsAsFactors = F, sep = "\t")
treatment_comp = read.csv(treatment_comp_file,  stringsAsFactors = F, sep = "\t")

d = merge(control_comp, treatment_comp,
          by = "id", suffixes = c("UNT","TRE"))
write.csv(d,output_file, quote = F, row.names = F)

colnames(d)
if (any(d$neg.p.valueUNT == 0.0)) {
  d[d$neg.p.valueUNT == 0.0,]$neg.p.valueUNT = min(d[d$neg.p.valueUNT > 0,]$neg.p.valueUNT) *0.9
}

d$UNT = -log10(d$neg.p.valueUNT)

if (any(d$neg.p.valueTRE == 0.0)) {
  d[d$neg.p.valueTRE == 0.0,]$neg.p.valueTRE = min(d[d$neg.p.valueTRE > 0,]$neg.p.valueTRE) *0.9
}
d$TRE = -log10(d$neg.p.valueTRE)
#plot(d$UNT, d$TRE)
if (all_genes) { 
start = max(d$TRE)
out = 0
count = 0 
while (out < number_of_genes_displayed & count < 100) {
  start = start  - 0.1
  out = dim(d[(d$UNT < -log10(0.01) & d$TRE > start),])[1]
  print(out)
  count = count + 1
}

smaller_set = d[(d$UNT < -log10(0.01) & d$TRE > start) |
                  d$id %in% gene_set,]
} else {
  smaller_set = d[d$id %in% gene_set,]
}

ggplot(d,aes(UNT,TRE, label = id)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$id %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("MAGeCK negative", norm_method, sep = " " )) +
  xlab(xlab_name) +
  ylab(ylab_name)
ggsave(negative_selection_plot, width = 3.5, height = 6)


if (any(d$pos.p.valueUNT == 0.0)) {
  d[d$pos.p.valueUNT == 0.0,]$pos.p.valueUNT = min(d[d$pos.p.valueUNT > 0,]$pos.p.valueUNT) *0.9
}
d$UNT = -log10(d$pos.p.valueUNT)

if (any(d$pos.p.valueTRE == 0.0)) {
  d[d$pos.p.valueTRE == 0.0,]$pos.p.valueTRE = min(d[d$pos.p.valueTRE > 0,]$pos.p.valueTRE) *0.9
}

d$TRE = -log10(d$pos.p.valueTRE)
#plot(d$UNT, d$TRE)

if (all_genes) {
start = max(d$TRE)
out = 0
count = 0 
while (out < number_of_genes_displayed & count < 100) {
  start = start  - 0.1
  out = dim(d[(d$UNT < -log10(0.01) & d$TRE > start),])[1]
  print(out)
  count = count + 1
}

smaller_set = d[(d$UNT < -log10(0.01) & d$TRE > start) |
                  d$id %in% gene_set,]
} else {
  smaller_set = d[d$id %in% gene_set,]
}


ggplot(d,aes(UNT,TRE, label = id)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$id %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("MAGeCK positive", norm_method, sep = " " )) +
  xlab(xlab_name) +
  ylab(ylab_name)
ggsave(positive_selection_plot, width = 3.5, height = 6)
