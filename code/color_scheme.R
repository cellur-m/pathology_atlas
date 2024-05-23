

#### gene classification
color_palette.gene_class <- c("cancer enriched" = "#e41a1c",
                              "group enriched" = "#ff9d00",
                              "cancer enhanced" = "#984ea3",
                               "low cancer specificity" = "#686868",
                              "not detected" = "#bebebe")
gene_class.levels <- factor(c("cancer enriched", "group enriched", "cancer enhanced", "low cancer specificity", "not detected"))





#### FD datasets
mutate_datasets = data.frame(
  order=c(1:10),
  abbr=c("LIHC", "KIRC", "PAAD", "COAD", "READ", "LUAD", "LUSC", "BRCA", "OV", "GBM"),
  HPA_V2=c("TCGA-LIHC", "TCGA-KIRC", "TCGA-PAAD", "TCGA-COAD", "TCGA-READ", "TCGA-LUAD", "TCGA-LUSC", "TCGA-BRCA", "TCGA-OV", "TCGA-GBM"),
  NewData=c("LIHC-FD", "KIRC-FD", "PAAD-FD", "COAD-FD", "READ-FD", "LUAD-FD", "LUSC-FD", "BRCA-FD", "OV-FD", "GBM-FD"),
  tissue=c("Liver", "Kidney", "Pancreas", "Intestine", "Intestine", "Lung", "Lung", "Breast", "Ovary", "Brain"),
  color=c("#D1CBE5", "#F9A266", "#96C08E", "#1280C4", "#1280C4", "#6AA692", "#6AA692", "#F8BDD7", "#F8BDD7", "#FFDD00")
)




