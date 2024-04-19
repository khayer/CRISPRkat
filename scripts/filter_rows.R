args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

counts_file = args[1]
print(counts_file)
output_file = args[2]
filter_cutoff = as.numeric(args[3])
print(filter_cutoff)
sample_names = unlist(strsplit(args[4],","))
print(sample_names)
print(length(sample_names))
#counts_file = "~/Dropbox (Work)/files/Abby/GeCKO_Sep_lab_meeting/code/data/median_norm_csgrnas.count.txt"
#output_file = "~/Desktop/tmp.txt"
#filter_cutoff = 5
#sample_names = c("d0_unt",
#                 "d15_unt", 
#                 "d15_ld", 
#                 "d15_hd")

data = read.csv(counts_file, sep="\t")
head(data)
print(dim(data))
data = data[rowSums(data[, colnames(data) %in% sample_names]) > filter_cutoff,]
print(dim(data))
write.table(data,output_file, sep="\t",row.names=FALSE,na='',quote = FALSE)
