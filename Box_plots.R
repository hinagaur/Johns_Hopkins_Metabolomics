# Loading necessary libraries
library(tidyverse)
library(magrittr)
library(readxl)
library(preprocessCore)
library(limma)
library(janitor)
library(umap)
library(Rtsne)
library(factoextra)
library(ggVennDiagram)
library(pheatmap)
library(xlsx)
library(KEGGREST)
library(curl)
library(ggplot2)
library(fgsea)
library(extrafont)
library(broom)
library(patchwork)
library(dplyr)
library(data.table)

#final data
df <- read.csv("/Users/hinagaur/Desktop/vip_metaboanalyst_data_3_and_2_cluster/clustering_3_vip_all.csv", check.names = FALSE)



#loading chem data
chem_data <- read_excel("/Users/hinagaur/Documents/KAMAL LAB/DrANN/chemical_annotation.xlsx")
meta_data <- read_excel("/Users/hinagaur/Documents/KAMAL LAB/DrANN/meta_data.xlsx")

#merging
df_new <- merge(df, meta_data, by.x = "PARENT_SAMPLE_NAME", by.y = "PARENT_SAMPLE_NAME") %>% select(-c("CLIENT_IDENTIFIER","NEG", "POLAR",
                                                                                                       "POS EARLY","POS LATE","AEROALLERGEN_SENSITIZA_1",
                                                                                                       "AGE_IN_MONTHS","BOX_NUMBER","CLIENT_MATRIX",
                                                                                                       "CLIENT_SAMPLE_ID","CLIENT_SAMPLE_NUMBER","ETHINICITY", "GENDER",
                                                                                                       "INHALED_STEROID", "RACE","SAMPLE_AMOUNT", "SAMPLE_AMOUNT_UNITS" ,"SAMPLE_BOX_LOCATION",
                                                                                                       "SAMPLE_DESCRIPTION" , "VOL_EXTRACTED")) %>% select(c(PARENT_SAMPLE_NAME, Cluster, GROUP_ID, everything()))


#unnest chem_data

chem_data_new <- chem_data %>% 
  mutate(KEGG = strsplit(as.character(KEGG), ",")) %>% 
  unnest(KEGG)

# correcting words
df_new$GROUP_ID <- gsub('Non_Sensitized','Non Sensitized',df_new$GROUP_ID)


# loading vip table
vip_table_3_clustering <- read.csv("/Users/hinagaur/Desktop/vip_metaboanalyst_data_3_and_2_cluster/cluster_3_all_vip_results/plsda_vip_table.csv")
colnames(vip_table_3_clustering)[1] <- "CHEMICAL_NAME"

vip_3_clust <- vip_table_3_clustering[1:15, ]


#Extracting metabolites
metab<- paste0( vip_3_clust$CHEMICAL_NAME, collapse="#")

#splitting
metab <- unlist(strsplit(metab,"#"))
print(metab)

comparisons <- list(c("Cluster 1", "Cluster 2", "Cluster 3"))

# Loop through each metabolite
for (m in metab) {
  print(paste("Processing metabolite:", m))
  
  # Create a boxplot for the current metabolite
  p <- ggplot(df_new, aes(x = Cluster, y = !!sym(m), fill = Cluster)) +
    geom_boxplot() +
    labs(title = m, x = "", y = "Normalized Concentration") +
    theme_bw() +
    geom_jitter(height = 0, width = 0.3) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c('Cluster 1' = 'blue', 'Cluster 2' = 'red', 'Cluster 3'= 'green')) +
    theme(
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 15, face = "bold", family = "Helvetica"),
      legend.title = element_text(size = 15),
      title = element_text(size = 15),
      legend.text = element_text(size = 15)
    ) +
    stat_compare_means(method = "t.test", comparisons = comparisons)
  
  # Save the plot
  filename <- paste("~/desktop/boxplot/boxplot_", gsub("[[:space:]/]", "_", m), ".png", sep = "")
  print(paste("Saving plot as:", filename))
  ggsave(filename, plot = p, width = 8, height = 6, units = "in")
  
  # Optionally, assign the plot to a variable in the R environment
  assign(paste("p", gsub(" ", "_", m), sep = ""), p)
}


###########
p <- ggplot(df_new, aes(x = GROUP_ID, y = xylose, fill = GROUP_ID)) +
  geom_boxplot() +
  labs(title = "xylose", x = "", y = "Normalized Concentration") +
  theme_bw() +
  geom_jitter(height = 0, width = 0.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c('Non Sensitized' = 'blue' , 'Sensitized' = 'red')) +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold", family = "Helvetica"),
    legend.title = element_text(size = 15),
    title = element_text(size = 15),
    legend.text = element_text(size = 15)
  )

p_with_pvalues <- p + stat_compare_means(method = "t.test", comparisons = list(c("Non Sensitized", "Sensitized")))
ggsave("~/desktop/xylose.png", plot = p_with_pvalues, width = 8, height = 6, units = "in")

