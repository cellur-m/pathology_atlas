### gene classify by specificity -----
rm(list=ls())

source("toolbox_pathology.R")
EnvironmentSetup()

options(scipen = 999)

# - main ---
data.raw <- read.csv("MeanTPM_hpa.v2.tsv", sep='\t') %>%
  column_to_rownames("ensg_id")

data.raw.long <- melt(as.matrix(data.raw)) %>%
  rename_with(., ~c('ensg_id', 'cancer', 'ave_TPM'))
head(data.raw.long)

cancers <- unique(data.raw.long$cancer)

# - gene classification ---
data.class <- hpa_gene_classification_cancer(data.raw.long, 'ave_TPM', "cancer", "ensg_id", 4, round(length(cancers)/3), det_lim = 1) %>%
  mutate(across(where(is.character), gsub, pattern="\\.", replacement = "\\-"))
table(data.class$spec_category)
table(data.class$dist_category)
table(data.class$n_expressed)

write.table(data.class, file = "hpa_gene_classification.v2.tsv", sep = '\t',
            col.names = T, row.names = F, quote = F)


# - visualization ---
# -- 1. numbers of gene class ---
source("color_scheme.R")

plot <-
  data.class %>%
  group_by(spec_category) %>%
  summarise(category_num=as.numeric(table(spec_category))) %>%
  ggplot(., aes(x=factor(spec_category, levels=gene_class.levels), y=category_num, fill=spec_category, label=category_num)) +
  geom_bar(stat="identity", width=0.7,
           position=position_dodge(width=0.8)) +
  scale_fill_manual(values=color_palette.gene_class, name="",
                    breaks = gene_class.levels) +
  geom_text(size=3, vjust=-0.2) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(angle=45, vjust=1, hjus=1),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,13000), expand = c(0, 0)) +
  xlab("") +
  ylab("Gene numbers") +
  ggtitle("Cancer gene specificity")







#############################################################
# -- 2. gene number in each cancer -----
# - functions --
statistic_by_category <-
  function(data, category){
    plot_data <- data %>%
      filter(spec_category==category)
    
    if (is.element("TRUE", grepl("\\;", plot_data$enriched_cancers))){
      # genes are classified into multiple cancers
      plot_data <- separate_rows(plot_data, enriched_cancers, sep=";")
    }
    
    plot.cancertype <- plot_data %>%
      group_by(enriched_cancers) %>%
      summarise(cancer_num=as.numeric(table(enriched_cancers))) %>%
      ggplot(., aes(x=reorder(enriched_cancers, -cancer_num), y=cancer_num, fill=enriched_cancers)) +
      geom_col()+
      scale_fill_manual(values = color_palette) +
      theme_bw()+
      theme(axis.text.x = element_text(angle=45, hjus=1, vjust=1),
            axis.title.y = element_text(size=8),
            panel.grid = element_blank(),
            legend.position = "none") +
      ylab(paste0(stringr::str_to_sentence(category), " gene number")) +
      xlab("")
    
    plot.cancertype
    
    ggsave(plot.cancertype, filename=paste0(category, " gene number.pdf"),
           width=10, height = 8, units = "cm", device = "pdf")
    
    plot.cancertype.label <- plot.cancertype +
      geom_text(aes(label=cancer_num), size=2, vjust=-0.2)
    
    ggsave(plot.cancertype.label, filename=paste0(category, " gene number_labeled.pdf"),
           width=10, height = 8, units = "cm", device = "pdf")
  }
# ---

# - main --
data.class <- read.csv('/hpa_gene_classification.tsv', sep='\t')

HPA_color <- read.csv("./colors_regions_pathology.tsv", sep='\t')
color_palette <- HPA_color %>%
  select(lable=dataset_name, color) %>%
  deframe()

for (cate in unique(data.class$spec_category)){
  # no cancer info in "low cancer specificity" and "not detected"
  message(cate)
  
  statistic_by_category(data.class, cate)
}
# ---
#############################################################
