##Loading necessary libraries
library(readxl)
library(dplyr)
library(data.table)
library(janitor)
library(reshape2)
library(fibroEset)
library(Biobase)
library(pheatmap)
library(ggrepel)
library(limma)
library(bestNormalize)
library(preprocessCore)


##setting maximum overlaps to infinite 
options(ggrepel.max.overlaps = Inf)

# LIMMA FOR M477 RENAL

# Setting working directory
setwd("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/")

## Metadata M477 renal metadata



# Reading Sample group info file
sample_group_renal <- read.delim("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/M477_VT/sample mapping files/rbc_mapping_c18neg_edit11152021_mauc.txt") %>% select(-c("Sample.ID",
                                                                                                                                                                                 "Batch","mauc_mes_bef","mauc_mes_dur",
                                                                                                                                                                                 "mauc_mes_aft", "mauc_mes_diff", "mauc_renal_bef", "mes_response_group",
                                                                                                                                                                                 "mauc_renal_dur", "mauc_renal_aft", "mauc_renal_diff", "mes_response_group_baseline_corrected",
                                                                                                                                                                                 "mes_residual_delta_mauc", "renal_residual_delta_mauc", "renal_response_group_baseline_corrected" )) %>% na.omit()

# Changing a colname to a different colname
colnames(sample_group_renal)[which(names(sample_group_renal) == "File.Name")] <- "Sample"

# Replacing 0 to no and 1 to yes
sample_group_renal$renal_response_group <- gsub('1', 'Yes', sample_group_renal$renal_response_group)
sample_group_renal$renal_response_group <- gsub('0','No',sample_group_renal$renal_response_group)

# Replacing YW in sample column to VT
# sample_group_renal$Sample <- gsub('YW_','VT_',sample_group_renal$Sample)

# Removing duplicates values under column
sample_group_renal <- sample_group_renal[!duplicated(sample_group_renal$Sample), ]


# getting all the sample names
sample.names <- sample_group_renal$Sample

# setting the first column as rownames names
rownames(sample_group_renal)<- sample_group_renal$Sample

#deleting extra column
sample_group_renal$Sample <- NULL


# Matrix data
# Loading feature table
feature_table <- read.delim("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/M477_VT/Extracted Data/c18neg/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt") %>% select(-c("mz.min","mz.max",
                                                                                                                                                                                                        "NumPres.All.Samples", "NumPres.Biological.Samples",
                                                                                                                                                                                                        "median_CV", "PeakScore", "Qscore", "Max.Intensity",
                                                                                                                                                                                                        "NumPres.All.Samples" ,"NumPres.Biological.Samples"))
# Joining mz and time with underscore
feature_table$mztime <- paste0(feature_table$mz, "__", feature_table$time)

# Reordering
feature_data <- feature_table %>% select(c("mztime", everything())) %>% select(-c("mz", "time"))

# Loading chem annotation table

chem_ann <- read.csv("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/M477_VT/Extracted Data/c18neg/xmsannotator_c18neg.csv")

# Joining mz and time with underscore
chem_ann$mztime <- paste0(chem_ann$mz, "__", chem_ann$time) 
chem_data <- chem_ann %>%  select(-c("chemical_ID", "Module_RTclust",
                                     "mz", "time", "MatchCategory", "theoretical.mz", "delta_ppm","Formula", "MonoisotopicMass",
                                     "Adduct", "mean_int_vec" )) %>% select(c("mztime", "Name"))

# Merging two columns

Final_df <- merge(feature_data, chem_data,  by = "mztime") %>% unique() %>% select(c("Name", everything())) %>% select(-c("mztime"))


#Subsetting to only keep specific columns
Final_df_new <- Final_df[, c("Name", sample.names)]

# Grouping and calculating median of redundant chem names
new_Final_df <- Final_df_new %>%
  group_by(Name) %>%
  summarise(across(everything(), ~median(., na.rm = TRUE)))


#converting df to matrix
mat <- as.matrix(new_Final_df)

#Converting first column to rownames
mat_new <- mat[,-1]
rownames(mat_new) <- mat[,1]


# check normalization
mat_numeric = apply(mat_new, 2, as.numeric)
rownames(mat_numeric)  = rownames(mat_new)

# boxplot(mat_numeric)

# Log transformation (added one cause of zero values)
log_data <- log(mat_numeric+1)
# boxplot(log_data)


# best normalization
scale = function(x){
  mean_ = mean(x, na.rm = T)
  sd_ = sd(x, na.rm = T)
  scaled = (x-mean_)/sd_
  return(scaled)
}

# standardize
scaled = apply(log_data, 2, scale)
rownames(scaled) = rownames(log_data)
# boxplot(scaled)

#Quantile normalization
quantiled = normalize.quantiles(log_data)
rownames(quantiled) = rownames(log_data)
colnames(quantiled) = colnames(log_data)
# boxplot(quantiled)

# #Chem data
# chem <- as.data.frame(rownames(mat_new))
# chem$chem_name <- chem$`rownames(mat_new)`
# 
# #converting to matrix
# rownames(chem) <- chem$`rownames(mat_new)`
# 
# #removing extra column
# chem$`rownames(mat_new)` <- NULL

## For Limma
# Quantiled - Normalized Matrix
# Chem - Chemical Annotation/Metadata
# Sample_group_mes - Sample Metadata


#Limma

#CREATING DESIGN MATRIX

design <- model.matrix(~0+ sample_group_renal$renal_response_group)

## renaming column names
colnames(design) <- c("RENAL_NO_RESPONSE","RENAL_YES_RESPONSE")

#The lmFit function is used to fit the model to the data.
#The result of which is to estimate the expression level in each of the groups that we specified.

fit <- lmFit(scaled, design)

#In order to perform the differential analysis, we have to define the contrast that we are interested in.
contrasts <- makeContrasts(RENAL_YES_RESPONSE - RENAL_NO_RESPONSE, levels=design)


fit2 <- contrasts.fit(fit, contrasts)

# Apply the empirical Bayes’ step to get our differential expression statistics and p-values.

fit2 <- eBayes(fit2)
topTable(fit2, coef=1)

#viewing top 10 differentially expressed metabolites
topTable(fit2)

##Saving entire results to an object

full_results <- topTable(fit2, number=Inf)

#Investigating up and down regulated metabolites
sig <- subset(full_results, P.Value < 0.05)
Down <- subset(full_results, P.Value < 0.05 & logFC < -2)
Up <- subset(full_results, P.Value < 0.05 & logFC > 2)

#Saving results
write.csv(full_results, "LIMMA_RESULTS/m477_renal_limma_log_scaled.csv")

