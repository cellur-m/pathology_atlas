# sig pathways between alive and dead patient groups
source("toolbox_pathology.R")
EnvironmentSetup()
library(ggridges)


# load clinical data
data.filter <- read.csv("Input.tsv", sep='\t') %>%
    FilterInfile(.)

samples <- data.filter$cli_data %>%
  tibble::column_to_rownames("sample_id") %>%
  mutate(status=as.factor(status)) %>%
  select(-c(stage, survival_time))
rownames(samples) <- gsub("-", "", rownames(samples))


# load data
load("network.RData")
Pathway_Activity <- Pathway_Activity %>%
  as.data.frame() %>%
  select(starts_with("KEGG")) %>%
  merge(., samples, by=0)

table(Pathway_Activity$status)

# find statistical different pathways
# -- Kolmogorov-Smirnov Tests
pathways <- colnames(Pathway_Activity)[which(grepl("KEGG", colnames(Pathway_Activity)))]
sample.alive <- rownames(samples)[which(samples$status=="alive")]
sample.dead <- rownames(samples)[which(samples$status=="dead")]

res <- tibble(Pathway=as.character(),
              p.ks=as.numeric())
for (i in pathways){
  x <- Pathway_Activity %>%
    filter(Row.names %in% sample.alive) %>%
    select(any_of(i)) %>% pull()
  y <- Pathway_Activity %>%
    filter(Row.names %in% sample.dead) %>%
    select(any_of(i)) %>% pull()

  
  ks.res <- ks.test(x, y, "two-sided")
  p.tmp <- ks.res$p.value
  
  res <- add_row(res,
                 Pathway=i,
                 p.ks=p.tmp)
}

res$p.adj <- p.adjust(res$p.ks, method = "BH")


top.ks.pathway <- res %>%
  filter(p.ks<0.05) %>%
  arrange(p.adj, desc=F) %>%
  slice_head(n=10) %>%
  mutate(order=c(1:n()))

