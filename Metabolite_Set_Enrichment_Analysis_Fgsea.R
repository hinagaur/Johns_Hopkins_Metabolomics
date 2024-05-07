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
library(KEGGREST)
library(curl)
library(ggplot2)
library(fgsea)
library(extrafont)
library(broom)


# Setting working directory
setwd("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/")

hmdb_data <- readRDS("/Users/hinagaur/Downloads/SMP_data 1.rds")

# pathway_name <- hmdb_data[[1]][[2]][[1]]
# hmdb_id <- hmdb_data[[1]][[6]]


# Making a named list of pathways as names and HMDB ID as values
new_list <- list()  # Create an empty list

# for (i in seq_along(hmdb_data)) {
#   if (length(hmdb_data[[i]]) >= 2 && length(hmdb_data[[i]][[2]]) >= 1) {
#     pathway_name <- hmdb_data[[i]][[2]][[1]]
#     hmdb_id <- hmdb_data[[i]][[6]]
#     new_list[[pathway_name]] <- hmdb_id
#   } else {
#     warning(paste("Invalid structure at index", i, "in hmdb_data. Skipping..."))
#   }
# }

for (i in seq_along(hmdb_data)) {
  if (length(hmdb_data[[i]]) >= 2 && length(hmdb_data[[i]][[2]]) >= 1) {
    pathway_name <- hmdb_data[[i]][[2]][[1]]
    hmdb_id <- hmdb_data[[i]][[6]]
    
    # Remove empty values from hmdb_id list
    hmdb_id <- hmdb_id[hmdb_id != ""]
    
    if (length(hmdb_id) > 0) {
      new_list[[pathway_name]] <- hmdb_id
    } else {
      warning(paste("Empty hmdb_id at index", i, "in hmdb_data. Skipping..."))
    }
  } else {
    warning(paste("Invalid structure at index", i, "in hmdb_data. Skipping..."))
  }
}

#Loading limma results for M407 renal

limma_res<- read.csv("LIMMA_RESULTS/m477_renal_limma_log_quantiled.csv")
# %>% subset( P.Value<0.05)

#Subsetting only for significant metabolites
# up <- subset(limma_res, P.Value<0.05 & logFC > 1)
# down <- subset(limma_res, P.Value<0.05 & logFC < -1)

#Changing colname of col contatining chem names
colnames(limma_res)[1] <- "Name"
# colnames(up)[1] <- "Name"
# colnames(down)[1] <- "Name"
#Loading chemical metadata
chem_data <- read.csv("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/M477_VT/Extracted Data/c18neg/xmsannotator_c18neg.csv")
chem <- chem_data[!duplicated(chem_data$chemical_ID), ]

#Adding two zeroes in the HMDB iDs
chem$chemical_ID <- str_replace(chem$chemical_ID, "B", "B00")

# Merging limma res with chem metadata
gseainput <- merge(chem, limma_res, by = "Name") %>% arrange(logFC)

# gseainput_up <- merge(chem, up, by = "Name") %>% arrange(logFC)
# gseainput_down <- merge(chem, down, by = "Name") %>% arrange(logFC)


# Making a named list
ranks_m477_renal <- setNames(gseainput$logFC, gseainput$chemical_ID)

# ranks_m407_renal_up <- setNames(gseainput_up$logFC, gseainput_up$chemical_ID)
# ranks_m407_renal_down <- setNames(gseainput_down$logFC, gseainput_down$chemical_ID)

#Fgsea

fgseaRes <- fgsea(pathways= new_list, stats=ranks_m477_renal, nperm=1000) %>% arrange(pval)

top_pathways <- subset(fgseaRes, pval <= 0.05) %>% as_tibble()
saveRDS(top_pathways, "/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/msea_results/msea_m477_renal_yw.rds")

top_pathways$pathway <- gsub("(?<=De Novo Triacylglycerol Biosynthesis).*", "", top_pathways$pathway, perl = TRUE)
top_pathways$pathway <- gsub("(?<=Adrenal Hyperplasia Type [3,5]).*", "", top_pathways$pathway, perl = TRUE)
# top_pathways$pathway <- gsub("(?=Adrenal Hyperplasia Type [3,5]).*", "", top_pathways$pathway, perl = TRUE)

top_pathways = top_pathways[!duplicated(top_pathways$pathway),]

# Filter out the top 20 upregulated and downregulated pathways
top_up <- top_pathways %>%
  filter(NES > 0) %>%
  top_n(10, wt = NES)

top_down <- top_pathways %>%
  filter(NES < 0) %>%
  top_n(10, wt = abs(NES)) %>% 
  slice(1:10)

top_10 <- bind_rows(top_up, top_down)
top_10$pathway <- gsub(".*(?=cblG Complementation Type)", "",  top_10$pathway, perl = TRUE)



# Plot the top 20 upregulated and downregulated pathways
p2_2clus <- ggplot(top_10, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = ifelse(NES > 0, "up", "down"))) +
  ylim(-2, 2) +
  coord_flip() +
  labs(x = " ", y = "Normalized Enrichment Score", title = " ") +
  scale_fill_manual(values = c("up" = "red", "down" = "blue"), name = "Direction") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 17, face = "bold", family = "Helvetica"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

# Print the plot
print(p2_2clus)


# p2_2clus = ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill = ifelse(NES > 0, "up", "down")))+
#   ylim(-2,2)+
#   coord_flip() +
#   labs( x = " ", y = "Normalized Enrichment Score", title = " ") +
#   scale_fill_manual(values = c("up" = "red", "down" = "blue"),
#                     name = "Direction") +
#   theme_minimal() +
#   theme(axis.title.x = element_text(size = 15, face= "bold"),
#         axis.text.x = element_text(size = 15, face ="bold"),     # Adjust the font size for x and y axis values
#         axis.text.y = element_text(size = 17, face= "bold", family = "Helvetica"),
#         legend.title = element_text(size = 15),    # Adjust the font size for the legend title
#         legend.text = element_text(size = 15))
# 
# #saving the msea results
# saveRDS(top_pathways, "/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/msea_results/msea_m477_renal_vt.rds")

