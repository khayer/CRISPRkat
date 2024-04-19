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
#gene_set = c("ATR","ATM","CHEK1","CHEK2")

# for testing only!!
#control_comp_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/pbnpa_total_day0_unt_day15_unt.csv"
#treatment_comp_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/pbnpa_total_day0_unt_day15_low_dox.csv"
#output_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/test_combined.csv"
#norm_method = "total"
#negative_selection_plot = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/negative.pdf"
#positive_selection_plot ="/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/positve.pdf"

ylab_name = gsub(".csv",x = tail(strsplit(treatment_comp_file,"/")[[1]], n=1),replacement = "")
ylab_name = gsub("^pbnpa_[A-Za-z]*_", ylab_name, replacement = "")

xlab_name = gsub(".csv",x = tail(strsplit(control_comp_file,"/")[[1]], n=1),replacement = "")
xlab_name = gsub("^pbnpa_[A-Za-z]*_", xlab_name, replacement = "")

control_comp = read.csv(control_comp_file,  stringsAsFactors = F)
treatment_comp = read.csv(treatment_comp_file,  stringsAsFactors = F)

d = merge(control_comp, treatment_comp,
          by = "Gene", suffixes = c("UNT","TRE"))
write.csv(d,output_file, quote = F, row.names = F)

colnames(d)
if (any(d$neg.pvalueUNT == 0.0)) {
  d[d$neg.pvalueUNT == 0.0,]$neg.pvalueUNT = min(d[d$neg.pvalueUNT > 0,]$neg.pvalueUNT) *0.9
}

d$UNT = -log10(d$neg.pvalueUNT)

if (any(d$neg.pvalueTRE == 0.0)) {
  d[d$neg.pvalueTRE == 0.0,]$neg.pvalueTRE = min(d[d$neg.pvalueTRE > 0,]$neg.pvalueTRE) *0.9
}
d$TRE = -log10(d$neg.pvalueTRE)
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
                  d$Gene %in% gene_set,]
} else {
  smaller_set = d[d$Gene %in% gene_set,]
}

ggplot(d,aes(UNT,TRE, label = Gene)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$Gene %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("PBNPA negative", norm_method, sep = " " )) +
  xlab(xlab_name) +
  ylab(ylab_name)
ggsave(negative_selection_plot, width = 3.5, height = 6)


if (any(d$pos.pvalueUNT == 0.0)) {
  d[d$pos.pvalueUNT == 0.0,]$pos.pvalueUNT = min(d[d$pos.pvalueUNT > 0,]$pos.pvalueUNT) *0.9
}
d$UNT = -log10(d$pos.pvalueUNT)

if (any(d$pos.pvalueTRE == 0.0)) {
  d[d$pos.pvalueTRE == 0.0,]$pos.pvalueTRE = min(d[d$pos.pvalueTRE > 0,]$pos.pvalueTRE) *0.9
}

d$TRE = -log10(d$pos.pvalueTRE)
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
                  d$Gene %in% gene_set,]
} else {
  smaller_set = d[d$Gene %in% gene_set,]
}


ggplot(d,aes(UNT,TRE, label = Gene)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$Gene %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("PBNPA positive ", norm_method, sep = " " )) +
  xlab(xlab_name) +
  ylab(ylab_name)
ggsave(positive_selection_plot, width = 3.5, height = 6)
