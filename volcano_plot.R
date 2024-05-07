

setwd("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/")

limma_res <- read.csv("limma_results_confounders/m407_mes_limma_bo_gw.csv") %>% select(-c("X"))
sig_test <- subset(limma_res, adj.P.Val < 0.05)
meta_chem <- read.csv("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/Metabolomics Datasets/M407_extracted_data_Joshua/c18neg/xmsannotator_c18neg.csv")

df <- merge(limma_res, meta_chem, by.x = "ID", by.y = "mz") %>% select(-c("chemical_ID","Annotation.confidence.score",
                                                                          "Annotation.raw.score", "Module_RTclust",
                                                                          "time", "MatchCategory", "theoretical.mz",
                                                                          "delta_ppm", "Formula", "MonoisotopicMass",
                                                                          "Adduct", "mean_int_vec"))

de <- df[!duplicated(df), ]

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$adj.P.Val< 0.1 & de$logFC > 1] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$adj.P.Val< 0.1 & de$logFC < -1] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$Name[de$diffexpressed != "NO"]

de <- read.csv("test.csv")
library(ggplot2)
library(ggrepel)

ggplot(data=de, aes(x=logFC, y=-log10(P.Value), col=diffexpressed, label=name2)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3) +  # Adjust the size value as needed
  scale_color_manual(values=c("blue", "black", "red")) +
  guides(color = guide_legend(title = NULL))


test <- de[order(de$adj.P.Val),]
test$name2 <- test$Name

write.csv(test,"test.csv")
# Specify the indices of cells to retain
indices_to_retain <- c(6734,10428, 10429,10430,10431, 13764, 13765, 13766, 13767,13768 )  # Example indices to retain

test$name2[!seq_along(test$name2) %in% indices_to_retain] <- NA

# Print the modified dataframe
print(df)

