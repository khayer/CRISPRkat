library(ScreenBEAM)

#Rscript {workflow.basedir}/scripts/run_screen_beam.R {input} {output} {params.name1} {params.name2} {params.sam_names1} {params.sam_names2} 2> {log}

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

r<-ScreenBEAM(

  ###input format
  input.file= counts_file #tab-separted file
  ,
  control.samples=sample_names1#column names of control samples
  ,
  case.samples=sample_names2#column names of case/treated samples
  ,
  control.groupname=g1_name #name your control group
  ,
  case.groupname=g2_name #name your case group
  ,

  ###data pre-processing
  data.type='NGS'#data type
  ,
  do.normalization=TRUE
  ,
  filterLowCount=TRUE
  ,
  filterBy = 'control'
  ,
  count.cutoff=4
  ,

  ###Bayesian computing
  nitt=15000,#number of MCMC iterations, use small number here for testing, please use larger number in real data, 15000 is default
  burnin=5000#number of burnin in MCMC sampling, 5000 is default

)
    
head(r)
###save your results
write.csv(r,file=file.path(output_file),row.names=FALSE,na='')
