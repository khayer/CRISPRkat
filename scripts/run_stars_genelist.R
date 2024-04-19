library(ggplot2)
library(ggrepel)

# Rscript {workflow.basedir}/scripts/run_screenbeam_genelist.R {input[0]} {input[1]} {output[0]} {params.norm_method} {output[1]} {output[2]} 

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

number_of_genes_displayed = 15

control_comp_file_1 = args[1]
control_comp_file_2 = args[2]
treatment_comp_file_1 = args[3]
treatment_comp_file_2 = args[4]
output_file = args[5]

norm_method = args[6]
negative_selection_plot = args[7]
positive_selection_plot = args[8]
essential = read.csv(args[9],sep="\t", header = F)
all_genes = args[10]

if (all_genes == "True") {
  all_genes = T
} else {
  F
}
gene_set = unlist(strsplit(args[11],","))
# for testing only!!
#control_comp_file_1 = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/06stats/day0_unt_vs_day15_unt_Score_STARSOutput_P.txt"
#control_comp_file_2 = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/06stats/day0_unt_vs_day15_unt_Score_STARSOutput_N.txt"
#treatment_comp_file_1 = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/06stats/day0_unt_vs_day15_low_dox_Score_STARSOutput_P.txt"
#treatment_comp_file_2 = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/06stats/day0_unt_vs_day15_low_dox_Score_STARSOutput_N.txt"
#output_file = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/STARS_test_combined.csv"
#norm_method = "control"
#negative_selection_plot = "/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/STARS_negative.pdf"
#positive_selection_plot ="/Users/hayerk/Documents/Matt/Abby/GeCKO/CRISPRkat_DEC/test/STARS_positve.pdf"
#essential = read.csv("~/Documents/Matt/Abby/GeCKO/essential_genes/core-essential-genes-sym_HGNCID",sep="\t", header = F)


test_if_skip_needed <- function(file_name) {
  n = 0
  # NEED TO TEST IF FIRST LINE STARTS WITH GENE
  first.line = readLines(file_name, n=1)
  if (startsWith(first.line,"p-value") ) {
    n = 1
  } 
  
  return(n)
}

ylab_name = gsub("_Score_STARSOutput_P.txt",x = tail(strsplit(treatment_comp_file_1,"/")[[1]], n=1),replacement = "")

xlab_name = gsub("_Score_STARSOutput_P.txt",x = tail(strsplit(control_comp_file_1,"/")[[1]], n=1),replacement = "")

 

control_comp_1 = read.csv(control_comp_file_1,  stringsAsFactors = F, sep = "\t", skip = test_if_skip_needed(control_comp_file_1))
control_comp_2 = read.csv(control_comp_file_2,  stringsAsFactors = F, sep = "\t", skip = test_if_skip_needed(control_comp_file_2))


control_comp_1$neg.p.value = control_comp_1$p.value
control_comp_1$pos.p.value = 1
control_comp_2$pos.p.value = control_comp_2$p.value
control_comp_2$neg.p.value = 1

control_comp = rbind(control_comp_1,control_comp_2)

treatment_comp_1 = read.csv(treatment_comp_file_1,  stringsAsFactors = F, sep = "\t", skip = test_if_skip_needed(treatment_comp_file_1))
treatment_comp_2 = read.csv(treatment_comp_file_2,  stringsAsFactors = F, sep = "\t", skip = test_if_skip_needed(treatment_comp_file_2))

treatment_comp_1$neg.p.value = treatment_comp_1$p.value
treatment_comp_1$pos.p.value = 1
treatment_comp_2$pos.p.value = treatment_comp_2$p.value
treatment_comp_2$neg.p.value = 1

treatment_comp = rbind(treatment_comp_1,treatment_comp_2)


d = merge(control_comp, treatment_comp,
          by = "Gene.Symbol", suffixes = c("UNT","TRE"))
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
                  d$Gene.Symbol %in% gene_set,]
} else {
  smaller_set = d[d$Gene.Symbol %in% gene_set,]
}

ggplot(d,aes(UNT,TRE, label = Gene.Symbol)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$Gene.Symbol %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("STARS negative", norm_method, sep = " " )) +
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
                  d$Gene.Symbol %in% gene_set,]
} else {
  smaller_set = d[d$Gene.Symbol %in% gene_set,] 
}


ggplot(d,aes(UNT,TRE, label = Gene.Symbol)) +
  geom_point(alpha = 0.5) +
  #geom_point(data=d[(d$UNT > 3 | d$TRE > 3) ,], col = "orange") +
  geom_point(data = d[d$Gene.Symbol %in% essential$V1,], aes(UNT,TRE), col = "blue", alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, col="lightgreen")  +
  geom_label_repel(data=smaller_set,fill = "white", col = "red", fontface = "bold",box.padding = 0.25,
                   point.padding = 0.3, alpha =0.7, size = 3 ) + theme_bw(base_size = 12) + 
  geom_point(data = smaller_set, col = "red", alpha = 0.6)+
  ggtitle(paste("STARS positive", norm_method, sep = " " )) +
  xlab(xlab_name) +
  ylab(ylab_name)
ggsave(positive_selection_plot, width = 3.5, height = 6)
