library("PBNPA")


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

counts_file = args[1]
print(counts_file)
output_file = args[2]
g1_name = args[3]
g2_name = args[4]
sample_names1 = unlist(strsplit(args[5],","))
sample_names2 = unlist(strsplit(args[6],","))



counts = read.csv(counts_file, sep ="\t", stringsAsFactors = F)
head(counts)

#counts = counts[counts$d0unt > 100,]
# compare d0unt to d15unt


indices_initial =  colnames(counts) %in% sample_names1
indices_final =  colnames(counts) %in% sample_names2
datlist = list()

if (length(sample_names2 ) > 1 ) {
  i = 1 
  for (k in sample_names2) {
    if (length(sample_names1) == 1) {
       r = 1
    } else {
      r = i
    }
    datlist[[i]] = data.frame("sgRNA" = counts$sgRNA,
                              "Gene" = counts$Gene,
                              "initial.count" = counts[, sample_names1[r]],
                              "final.count" = counts[, k])
    i = i +1
  }
  
} else {
  datlist = list(data.frame("sgRNA" = counts$sgRNA,
                            "Gene" = counts$Gene,
                            "initial.count" = counts[, indices_initial],
                            "final.count" = counts[, indices_final]))


}
#head(datlist)



result = PBNPA(datlist)

write.csv(result$final.result, output_file, row.names=FALSE,na='')
