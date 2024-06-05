EnvironmentSetup <- function(Update){
  if (!hasArg("Update")) {
    Update = F
  }
  if (Update == T) {
    install.packages("purrr")
    install.packages("tibble")
    install.packages("ggpubr")
    install.packages("ggplot2")
    install.packages("ggthemes")
    install.packages("dplyr")
    install.packages("tidyr")
    install.packages("tidyverse")
    install.packages("venn")
    install.packages("magrittr")
    install.packages("imputeTS")
    install.packages("xlsx")
    install.packages("reshape2")
    
    # source("https://bioconductor.org/biocLite.R")
    # biocLite("DESeq2")
    # biocLite("limma")
    # biocLite("DESeq2")
    # #biocLite("BiocUpgrade")
    # biocLite("piano", dependencies=TRUE)
    # biocLite("GO.db")
    # biocLite("Biobase")
    # biocLite("biomaRt")
  }
  
  library("purrr")
  library("tibble")
  library("ggpubr")
  library("ggplot2")
  library("ggthemes")
  library("dplyr")
  library("tidyr")
  library("tidyverse")
  library("venn")
  library("magrittr")
  library("imputeTS")
  library("readxl")
  library("xlsx")
  library("reshape2")
  library("ComplexHeatmap")
  library("corrplot")
  
  select <- dplyr::select
  
}



#####
# format: first col: sample ID;   colnames:gene ensemble ID
FilterInfile <- function(rawdata){
  # -- extract EXP file
  exp_index <- min(grep('ENSG', colnames(rawdata)))
  exp_data <- rawdata[, exp_index:ncol(rawdata)] %>%
    magrittr::set_rownames(rawdata[, 1])
  ## -- delete donors with all NA tpm
  ind <- apply(exp_data, 1, function(x) all(is.na(x)))
  exp_data <- exp_data[!ind, ]
  message(paste("Found", length(rownames(exp_data)[ind]), "sample(s) with all NA TPM values"))
  ## -- replace NA value with 0
  exp_data[is.na(exp_data)] <- 0
  ## -- keep TPM > 1
  sample_keep <- rownames(exp_data)
  exp_data<- apply(exp_data, 2, as.numeric)
  rownames(exp_data) <- sample_keep
  exp_keep <- as.data.frame(exp_data[, which(colMeans(exp_data)>1)])
  
  # -- extract clinical information
  cli_data <- rawdata[, 1:exp_index-1]
  cli_data <- cli_data[, c("sample_id", "status", "survival_time")]

  result <- list()
  result$exp_keep <- exp_keep
  result$cli_data <- cli_data
  return(result)
}
#####



#####
# given two database names and their exp profiles
# format: rownames:sample;  colnames:gene ensemble ID
MeanEXPCompare <- function(name1, exp_keep1, name2, exp_keep2){
  # -- calculate gene mean expression
  exp_keep1['mean_exp_1', ] <- colMeans(exp_keep1)
  exp_keep2['mean_exp_2', ] <- colMeans(exp_keep2)
  # -- get log2(mean exp) table
  mutate_genes <- intersect(colnames(exp_keep1), colnames(exp_keep2))
  diff_exp <- as.data.frame(t(rbind(exp_keep1['mean_exp_1', mutate_genes],
                                    exp_keep2['mean_exp_2', mutate_genes]))) %>%
    mutate(mean_exp_1 = log2(mean_exp_1),
           mean_exp_2 = log2(mean_exp_2))
  # -- plot with spearman correlation
  p = ggscatter(diff_exp, x='mean_exp_1', y='mean_exp_2', size=0.1, color = "#00203F",
                add.params = list(color = "black", fill = "gray20", size=0.1),
                add = 'reg.line',  conf.int = TRUE, conf.int.level = 0.95,
                xlim=c(0, 15), ylim=c(0, 15),
                xlab = paste0(name1, ' log2(TPM)'), ylab = paste0(name2, ' log2(TPM)')) +
    stat_cor(method='spearman', size =  4,  r.accuracy = 0.01,
             label.x = 1, label.y = 14) +
    geom_hline(yintercept = 0, color = 'gray', linetype='dashed', size = 0.5) +
    geom_vline(xintercept = 0, color = 'gray', linetype='dashed', size = 0.5) +
    theme(axis.text=element_text(size=1),
          panel.border = element_rect(colour='black', size=1, fill=NA)) +
    ggtitle("Mean expression level correlation")
  print(p)
  return(p)
}
#####


#####
# compare clinical differences
CliInfoCompare <- function(name1, cli_data1, name2, cli_data2){
  cli_data1[, 'cohort'] <- name1
  cli_data2[, 'cohort'] <- name2
  cli_merge <- rbind(cli_data1, cli_data2) %>%
    mutate(status=gsub('dead', 'deceased', status))
  
  p = ggplot(cli_merge, aes(x=status, y=survival_time, fill=cohort)) +
    scale_fill_manual(values=c("black", "black")) +
    geom_boxplot(alpha=0) +
    geom_point(aes(color=cohort), size=0.2, shape =20,position=position_jitterdodge()) +
    scale_color_manual(values=c("#B02179", "#07A082")) +
    stat_compare_means(method='t.test', size=2.8, aes(label = ..p.signif..)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = 'black', size=1, fill=NA),
          legend.background=element_blank()) +
    ggtitle("Survival information comparison") +
    xlab("Survival status") + ylab("Survival time / day")
  print(p)
  return(p)
}
#####


#####
# compare clinical differences
CliInfoCompare_violin <- function(name1, cli_data1, name2, cli_data2){
  cli_data1[, 'cohort'] <- name1
  cli_data2[, 'cohort'] <- name2
  cli_merge <- rbind(cli_data1, cli_data2) %>%
    mutate(status=gsub('dead', 'deceased', status))

  p = ggplot(cli_merge, aes(x=status, y=survival_time, fill=cohort)) +
    scale_fill_manual(values=c("black", "black")) +
    geom_violin(alpha=0) +
    geom_point(aes(color=cohort), size=0.2, shape =20,position=position_jitterdodge()) +
    scale_color_manual(values=c("#B02179", "#07A082")) +
    stat_compare_means(method='t.test', size=2.8, aes(label = ..p.signif..)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = 'black', size=1, fill=NA),
          legend.background=element_blank()) +
    ggtitle("Survival information comparison") +
    xlab("Survival status") + ylab("Survival time / day")
  print(p)
  return(p)
}
#####



#####
# compare KM genes with same direction
PlotKmCorr <- function(data){
  p = ggscatter(data, x = 'coef_1', y = 'coef_2', size = 2, color='Type',
                palette= c('gray', '#08326E'),
                add.params = list(color = "blue", fill = "gray20", size=0.7),
                add = 'reg.line',  conf.int = TRUE, conf.int.level = 0.95,
                xlim=c(-2, 2), ylim=c(-2, 2)) +
    stat_cor(method='spearman', size =  7,
             label.x = -2, label.y = 2) +
    geom_hline(yintercept = 0, color = 'gray', linetype='dashed', size = 0.5) +
    geom_vline(xintercept = 0, color = 'gray', linetype='dashed', size = 0.5) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_rect(colour = 'black', size=1, fill=NA),
          legend.background=element_blank(),
          legend.position=c(0.9, 0.1),
          legend.text = element_text(size=10)) +
    ggtitle("Kaplan-Meier coefficient correlation") +
    xlab("Survival status") + ylab("Survival time / day")
  print(p)
  return(p)
}
#####


LonglistToLongDataframe <- function(list){
  require(dplyr)
  require(magrittr)
  df_result <- data.frame(matrix(NA, ncol=2)) %>%
    set_colnames(c("name", "element"))
  
  for (tmp in names(list)){
    tmp_df <- data.frame(name=tmp,
                         element=list[tmp]) %>%
      set_colnames(c("name", "element"))
    df_result <- rbind(df_result, tmp_df)
    
    rm(tmp_df)
  }
  df_result <- df_result[complete.cases(df_result), ]
  message(dim(df_result)[1])
  return(df_result)
}



GenerateSurvInput <- function(gene,cancerType,data) {
  LivingDays = as.numeric(data$survival_time)
  SurvInput= as.data.frame(LivingDays)
  SurvInput$EXP = as.numeric(data[,gene])
  SurvInput$DeadInd = data$status %in% c('1')
  SurvInput$SurvObj <- with(SurvInput, Surv(LivingDays, DeadInd))
  return(SurvInput)
}




generateKMplot <- function(SurvInput,quantCut,outFile,figureTitle) {
    
  if (!hasArg("quantCut")) {
    cutoffs = sort(unique(SurvInput$EXP))
    Percentile20 = quantile(SurvInput$EXP,0.2, na.rm=T)
    Percentile80 = quantile(SurvInput$EXP,0.8, na.rm=T)
    cutoffs = cutoffs[cutoffs>Percentile20&cutoffs<=Percentile80]
    Pvalue = 1
    cutoff = 0
    mini = 0
    Pvalues = cutoffs
    for (i in 1:length(cutoffs)){
      tempdata = SurvInput
      tempdata$EXP[tempdata$EXP<cutoffs[i]]=0
      tempdata$EXP[tempdata$EXP>=cutoffs[i]]=1
      survOut = survfit(SurvObj ~ EXP,tempdata)
      v <- as.vector(tempdata)
      if (length(unique(v)) == 1){
        res = survdiff(SurvObj ~ EXP, SurvInput)
        logRankP = 1
        cutoff<-NA
        coef<-NA
        result = as.data.frame(logRankP)
        result$EXPcut = cutoff
        result$coef=coef
      } else {
        res.km2 <- survdiff(SurvObj~ EXP, tempdata, rho=0)
        icutoff = cutoffs[i]
        iPvalue = pchisq(res.km2$chisq,length(res.km2$n)-1,lower.tail = FALSE)
        Pvalues[i] = iPvalue
        if (iPvalue < Pvalue) {
          cutoff = icutoff
          Pvalue = iPvalue
          mini = i
      }
      }
      
    }
  } else {
    cutoff = quantile(SurvInput$EXP,quantCut)
  }
  SurvInput$EXP[SurvInput$EXP<cutoff]=0
  SurvInput$EXP[SurvInput$EXP>=cutoff]=1
  survOut = survfit(SurvObj ~ EXP,SurvInput)
  
  
  pdf(file = paste(outFile,"pdf",sep = "."),width=5.8,height=6)
  plot(survOut, col=c("red","black"), mark.time=T, cex=1.4,xlab="Time (year)",xscale = 365,lty =1, ylab = "Survival Probability",las=1, cex.lab=1.4)
  group1legend = paste("Low expression",paste("(n =",as.character(sum(SurvInput$EXP==0)),")",sep = ""), sep = " ")
  group2legend = paste("High expression",paste("(n =",as.character(sum(SurvInput$EXP==1)),")",sep = ""), sep = " ")
  if (!hasArg("figureTitle")) {
    #title("ASS1")
  } else {
    title(figureTitle)
  }
  legend(
    "bottomleft",
    legend=c(group1legend,group2legend),
    col=c("red","black"),
    horiz=FALSE,
    lty=1,
    bty='n',
    cex=1.4)

  res = survdiff(SurvObj ~ EXP, SurvInput)
  logRankP = 1 - pchisq(res$chisq, length(res$n)-1)
  legend("topright",legend =c(paste0("P=",as.character(format(logRankP,scientific = TRUE,digits = 3)))),text.font=2,bty="n",cex=1.4)
  dev.off()
  sum_result<-summary(coxph(SurvObj ~ EXP, SurvInput))
  coef<-sum_result$coefficients[1]
  result = as.data.frame(logRankP)
  result$EXPcut = cutoff
  result$coef=coef
  return(result)
}

