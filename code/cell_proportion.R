source("./toolbox_pathology.R")
EnvironmentSetup()

source("./Deconvolution_functions.R")


### set input and result path
dir.create("./results/")


### load bulk data
dataBulk <- read.csv(dataBulk.file, sep = '\t')
dataBulk <- dataBulk %>%
  select(sample_id, starts_with("ENSG")) %>%
  column_to_rownames("sample_id") %>%
  t() %>% as.data.frame() %>%
  rownames_to_column("ensg_id")


### load sc data
dataSC <- read.csv(dataSC.file, sep='\t', check.names=F) %>%
  column_to_rownames("cellType")
  #column_to_rownames("X")
labels <- sapply(strsplit(colnames(dataSC), "\\."), function(x) x[1])

# deconvolution
Signature <- buildSignatureMatrixMAST(dataSC, labels, "results")

allCounts_DWLS<-NULL
allCounts_OLS<-NULL
allCounts_SVR<-NULL
for(j in 1:(dim(dataBulk)[2]-1)){
  S<-Signature
  Bulk<-dataBulk[,j+1] # extract every sample
  names(Bulk)<-dataBulk[,1]
  Genes<-intersect(rownames(S),names(Bulk))
  B<-Bulk[Genes]
  S<-S[Genes,]
  solDWLS<-solveDampenedWLS(S,B)
  
  allCounts_DWLS<-cbind(allCounts_DWLS,solDWLS)

}

colnames(allCounts_DWLS) <- colnames(dataBulk)[2:dim(dataBulk)[2]]


# save solutions
save(allCounts_DWLS,file="results/allCounts_DWLS.RData")



