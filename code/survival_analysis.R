rm(list=ls())
source("toolbox_pathology.R")
EnvironmentSetup()


infile <- read.csv("example_SurvivalAnalysisInput.tsv", sep='\t')

columns <- colnames(infile)
first_gene_col <- 4
gene_list <- colnames(infile)[first_gene_col:ncol(infile)]

p_list<-NULL
cutoff_list<-NULL
coef_list<-NULL
for (j in 1:length(gene_list)){
  print(j)
  output<-GenerateSurvInput(gene_list[j],cancer_name,infile)
  setwd(path_out_pdf)
  result<-generateKMplot(output,outFile=paste0(cancer_name,"_",gene_list[j]))
  log_rank_p<-as.matrix(result$logRankP)
  cut_off<-as.matrix(result$EXPcut)
  coef<-as.matrix(result$coef)
  
  p_list<-rbind(p_list,log_rank_p)
  cutoff_list<-rbind(cutoff_list,cut_off)
  coef_list<-rbind(coef_list,coef)
  rm(log_rank_p,cut_off,coef)
}
result<-cbind(as.matrix(gene_list),coef_list,cutoff_list,p_list)
colnames(result)<-c("ensemble","coef","cutoff","p")
result <- as.data.frame(result)

