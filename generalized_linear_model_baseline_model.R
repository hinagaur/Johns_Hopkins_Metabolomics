## Seibi code

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
library(stringr)
library(openxlsx)
install.packages("readr")
library(readr)
library(car)

# Setting working directory
setwd("/Users/hinagaur/Documents/KAMAL LAB/untargeted_metabolomics/")

# transfusion date information
demo1 = read_excel("transfusion_date_dr_joshua/samplekey_179 transfused units sent 6-17-20.xlsx", skip = 2) %>%
  dplyr::rename(
    reci_patient_id = `Study ID`,
    donor_unit_id = `Donor Unit ID`,
    donor_sample_id = `sample ID`,
    trans_date = `Date of Transfusion`,
    irrad = `irrad/  nonirrad`,
    wk = `wk#`) %>%
  mutate(trans_date = as.numeric(trans_date)) %>%
  mutate(trans_date = as.Date(trans_date, origin = "1899-12-30")) %>%
  dplyr::select(reci_patient_id, donor_unit_id, donor_sample_id, trans_date, irrad, wk) %>%
  filter(!is.na(trans_date))


demo2 = read_excel("transfusion_date_dr_joshua/samplekey_331 samples-all_transfusins.xlsx", skip = 2)%>%
  dplyr::rename(
    reci_patient_id = id,
    donor_unit_id = `DonorUnitId031`,
    donor_sample_id = `sample ID`,
    trans_date = `Date of Transfusion`,
    irrad = `irrad/  nonirrad`,
    wk = `wk#`) %>%
  mutate(trans_date = as.numeric(trans_date)) %>%
  mutate(trans_date = as.Date(trans_date, origin = "1899-12-30")) %>%
  dplyr::select(reci_patient_id, donor_unit_id, donor_sample_id, trans_date, irrad, wk) %>%
  filter(!is.na(trans_date))

# split reci patient id
temp_list = str_split(demo2$reci_patient_id, "")
new_list = lapply(temp_list, function(x){
  paste0(x[1],"-", x[2],"-", x[3],x[4],x[5],x[6], "-", x[7])
})   
# add
demo2 = demo2 %>% mutate(reci_patient_id = do.call(rbind, new_list)[,1])


demo3 = read_excel("transfusion_date_dr_joshua/samplekey_list of 120 samples -units.xlsx", skip = 2)%>%
  dplyr::rename(
    reci_patient_id = `Study ID`,
    donor_unit_id = `Donor Unit ID`,
    donor_sample_id = `sample ID`,
    trans_date = `Date of Transfusion`,
    irrad = `irrad/  nonirrad`,
    wk = `wk#`) %>%
  mutate(trans_date = as.Date(trans_date)) %>%
  dplyr::select(reci_patient_id, donor_unit_id, donor_sample_id, trans_date, irrad, wk) %>%
  filter(!is.na(trans_date))


demo4 = read_excel("transfusion_date_dr_joshua/samplekey_Transfused Units (52) sent 6-29-20.xlsx", skip =2)%>%
  dplyr::rename(
    reci_patient_id = `Study ID`,
    donor_unit_id = `Donor Unit ID`,
    donor_sample_id = `sample ID`,
    trans_date = `Date of Transfusion`,
    irrad = `irrad/  nonirrad`,
    wk = `wk#`) %>%
  mutate(trans_date = as.Date(trans_date)) %>%
  dplyr::select(reci_patient_id, donor_unit_id, donor_sample_id, trans_date, irrad, wk) %>%
  filter(!is.na(trans_date))


trans_date_table = bind_rows(demo1, demo2, demo3, demo4) %>%
  mutate(donor_sample_unit_id = str_extract(donor_sample_id, regex("\\w*"))) %>%
  distinct(.)


## A tibble: 680 × 7
#   reci_patient_id donor_unit_id donor_sample_id trans_date irrad wk    donor_sample_unit_id
#   <chr>           <chr>         <chr>           <date>     <chr> <chr> <chr>               
# 1 1-2-0100-2      W200319465333 035RC21(1)      2019-02-26 irrad W-1   035RC21             
# 2 1-2-0100-2      W200319460693 035RC22(1)      2019-03-04 irrad W-1   035RC22             
# 3 1-2-0100-2      W204219030498 035RC23(1)      2019-03-11 irrad W-3   035RC23             
# 4 1-2-0100-2      W204219030498 035RC23(2)      2019-03-14 irrad NA    035RC23             
# 5 1-2-0100-2      W200319460392 035RC24(1)      2019-03-18 irrad W-4   035RC24             
# 6 1-2-0100-2      W200319482763 035RC25(1)      2019-04-01 irrad W-6   035RC25             
# 7 1-2-0100-2      W200319477961 035RC26(1)      2019-04-15 irrad W-8   035RC26             
# 8 1-2-0100-2      W200319477961 035RC26(2)      2019-04-15 irrad NA    035RC26   







#############################
# donor blood baseline
donors = read_csv("donor_id_information.csv") %>%
  dplyr::rename(
    reci_patient_id = id,
    donor_unit_id = DonorUnitId031,
    date_trans = DateTransfusion031)


donors = donors %>% arrange(donor_unit_id) %>%
  mutate(date_blood_draw = date_trans - age_of_blood) %>%
  mutate(date_irradiation = date_blood_draw + age_at_irradiation)  %>%
  dplyr::select(donor_unit_id, date_blood_draw, date_irradiation) %>%
  filter(!duplicated(donor_unit_id))

## A tibble: 184 × 3
#   donor_unit_id date_blood_draw date_irradiation
#   <chr>         <date>          <date>          
# 1 W115117212991 2017-09-09      2017-09-09      
# 2 W115117216593 2017-08-29      2017-09-07      
# 3 W115117241624 2017-09-03      2017-09-14      
# 4 W115119150896 2019-10-10      2019-10-18      
# 5 W115120018758 2020-02-06      2020-02-06      
# 6 W115120018758 2020-02-06      2020-02-06      
# 7 W115120045425 2020-02-05      2020-02-11      


trans_data = trans_date_table %>%
  left_join(donors, by = "donor_unit_id") %>%
  filter(!is.na(date_blood_draw))

dups = trans_data %>% filter(duplicated(donor_unit_id)) %>% pull(donor_unit_id)
trans_data %>% filter(donor_unit_id %in% dups) %>%
  arrange(donor_unit_id) %>%
  filter(!is.na(date_irradiation)) %>%
  print(n=100)

#   reci_patient_id donor_unit_id donor_sample_id trans_date
#   <chr>           <chr>         <chr>           <date>    
# 1 1-3-0102-1      W115117165487 022RC36(1)      2017-11-30
# 2 1-3-0102-1      W115117165487 022RC36         2017-11-30
# 3 1-3-0102-2      W115117165487 023RC35(1)      2017-11-26
# 4 1-3-0102-2      W115117165487 023RC35         2017-11-26

# same donor blood unit is transfused to different patients, in this case, donor_sample_id are different!!
# so. donor_unit_id, donor_sample_id, trans_date gets unique combination

# trans_data
# n = 419

# to merge the recipient demo/metabolomics, donor_unit_sample_id is sufficient
trans_data = trans_data %>% dplyr::select(donor_sample_unit_id, date_blood_draw, date_irradiation) %>%
  distinct(.)
# n = 185



#######################
# recipient demographics
reci_demo = read_excel("demographic_data_old.xlsx", sheet = "m407_mes") %>%
  dplyr::rename(
    blinded_id = Sample,
    reci_sample_id = Sample.ID,
    reci_patient_id = ID)
# split reci patient id
temp_list = str_split(reci_demo$reci_patient_id, "")
new_list = lapply(temp_list, function(x){
  paste0(x[1],"-", x[2],"-", x[3],x[4],x[5],x[6], "-", x[7])
})   
# add
reci_demo = reci_demo %>% mutate(reci_patient_id = do.call(rbind, new_list)[,1]) %>%
  # transfusion date
  mutate(trans_date = str_split(reci_sample_id, "_", simplify = TRUE)[,3]) %>%
  mutate(trans_date = as.POSIXct(trans_date, format = "%m%d%Y", tx = "EST"))



###################################
# outcome
response_meta = read.table("metab_annotationBySamples/m407/transfusion_mapping_hilicpos_mauc.txt", sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  # remove duplication based on all columns
  distinct(.) %>%
  dplyr::rename(
    donor_sample_id = File.Name,
    reci_sample_id = SampleID) %>%
  mutate(donor_sample_id = str_replace(donor_sample_id, regex("\\w{2}_"), "")) %>%
  mutate(reci_patient_id =  str_extract(reci_sample_id, regex("\\d{7}"))) %>%
  mutate(donor_sample_id = str_split(reci_sample_id, "_", simplify = TRUE)[,2])
# split reci patient id
temp_list = str_split(response_meta$reci_patient_id, "")
new_list = lapply(temp_list, function(x){
  paste0(x[1],"-", x[2],"-", x[3],x[4],x[5],x[6], "-", x[7])
})   
# add
response_meta = response_meta %>% mutate(reci_patient_id = do.call(rbind, new_list)[,1]) %>%
  mutate(trans_date = str_split(reci_sample_id, "_", simplify = TRUE)[,3]) %>%
  mutate(trans_date = as.POSIXct(trans_date, format = "%m%d%Y", tx = "EST"))




# outcomes in hilic pos and c18 neg should be same?
response_meta2 = read.table("metab_annotationBySamples/m407/transfusion_mapping_c18neg_mauc.txt", sep = "\t", header = TRUE) %>%
  as_tibble() %>%
  # remove duplication based on all columns
  distinct(.) %>%
  dplyr::rename(
    donor_sample_id = File.Name,
    reci_sample_id = Sample.ID) %>%
  mutate(donor_sample_id = str_replace(donor_sample_id, regex("\\w{2}_"), "")) %>%
  mutate(reci_patient_id =  str_extract(reci_sample_id, regex("\\d{7}"))) %>%
  mutate(donor_sample_id = str_split(reci_sample_id, "_", simplify = TRUE)[,2])
# split reci patient id
temp_list = str_split(response_meta2$reci_patient_id, "")
new_list = lapply(temp_list, function(x){
  paste0(x[1],"-", x[2],"-", x[3],x[4],x[5],x[6], "-", x[7])
})   
# add
response_meta2 = response_meta2 %>% mutate(reci_patient_id = do.call(rbind, new_list)[,1]) %>%
  mutate(trans_date = str_split(reci_sample_id, "_", simplify = TRUE)[,3]) %>%
  mutate(trans_date = as.POSIXct(trans_date, format = "%m%d%Y", tx = "EST"))



response_meta_temp = response_meta %>% select(reci_patient_id, trans_date, mes_response_group) %>%
  filter(!is.na(mes_response_group))

response_meta_temp2 = response_meta2 %>% select( reci_patient_id, trans_date, mes_response_group, mauc_mes_aft, mauc_mes_bef) %>%
  filter(!is.na(mes_response_group))

# response_meta_temp %>% full_join(response_meta_temp2, by = c("reci_patient_id", "trans_date")) %>%
# print(n=100)
# different. The values are the same, but some row are missing in either, and the other had a value
# so, it is better to take an advantage using both


# outcomes
# outcome = response_meta_temp %>% full_join(response_meta_temp2, by = c("reci_patient_id", "trans_date")) %>%
#   mutate(mes_response =
#            ifelse(is.na(mes_response_group.x),mes_response_group.y,
#                   ifelse(is.na(mes_response_group.y),mes_response_group.x, mes_response_group.x))) %>%
#   dplyr::select(reci_patient_id, trans_date, mes_response)

outcome = response_meta_temp2 %>%   dplyr::select(reci_patient_id, trans_date, mes_response_group, mauc_mes_aft, mauc_mes_bef) 


# merge demo
reci_demo = reci_demo %>% full_join(outcome, by = c("reci_patient_id", "trans_date"))

reci_demo %>%  full_join(outcome, by = c("reci_patient_id", "trans_date"))  %>%
  filter(is.na(DOB)) %>%
  print(n=100, width = Inf)

#  blinded_id reci_sample_id reci_patient_id DOB    sex   birthweight race  gestage_weeks trans_date         
#  <chr>      <chr>          <chr>           <dttm> <chr>       <dbl> <chr>         <dbl> <dttm>             
#1 NA         NA             1-2-0002-1      NA     NA             NA NA               NA 2017-05-03 00:00:00
#2 NA         NA             1-2-0100-1      NA     NA             NA NA               NA 2019-03-11 00:00:00
#3 NA         NA             1-3-0309-2      NA     NA             NA NA               NA 2019-03-02 00:00:00
#4 NA         NA             1-3-0329-1      NA     NA             NA NA               NA 2019-04-20 00:00:00
#  mes_response
#         <int>
#1            1
#2            0
#3            0
#4            1
# for patients are missing for information.





# to get blood baseline information, donor_unit_sample_id is sufficient.
reci_demo = reci_demo %>%
  mutate(donor_sample_id = str_split(reci_sample_id, "_", simplify = TRUE)[,2]) %>%
  mutate(donor_sample_unit_id = str_extract(donor_sample_id, regex("\\w*"))) %>%
  mutate(trans_date =as.Date(trans_date)) %>%
  left_join(trans_data, by = "donor_sample_unit_id")

# write.xlsx(reci_demo, "reci_demo_c18neg.xlsx")
# saveRDS(reci_demo, "recipient_demo_c18_neg.rds")




############extraaa just checking the information for 

# heli <- subset(reci_demo, is.na(blinded_id)) %>%  select(c("reci_patient_id", "trans_date", "mes_response"))

# 
# heli$reci_patient_id <- gsub("-","", heli$reci_patient_id) # Patient info present!

######### merging info to get donor blood info

# Loading transfusion dates from files that dr joshua had sent
trans_info <- read.csv("transfusion_date_dr_joshua/samplekey_331 samples-all_transfusins_copy_1.csv")
corr_demo <- read.csv("corrected_demo_1.csv")
trans_info = trans_info[-1,]

# Making first row as header 
names(trans_info) <- trans_info[1,]

#Deleting first row
trans_info = trans_info[-1,]

# splitting Sample.ID column 
corr_demo[c('id', 'sample ID' , 'date_trans', 'type', 'number')] <- str_split_fixed(corr_demo$Sample.ID, '_', 5)
colnames(corr_demo)[4] <- "DateTransfusion031"

#merging column based on three columns id, sampleID, DateTransfusion031
merged_obj <- merge(corr_demo, trans_info, by = c("id", "sample ID", "DateTransfusion031"))
# write.xlsx(merged_obj,"donor_blood_info_req.xlsx")

# Checking duplicated blood IDS
duplicated(merged_obj$DonorUnitId031)

# final merging
final_merge <- merge(merged_obj, reci_demo, by.x = "Sample.ID", by.y = "reci_sample_id") %>% unique() 

# There are duplicates in terms of the blood draw date and date of irradiation, so we remove the later blood draw and irradiation

Final <- final_merge %>% 
  arrange(date_blood_draw,date_irradiation) %>% 
  filter(!duplicated(Sample.ID)) 


# RUN GLM ON ALL MEATBOLITES

# FIRST EXTRACT ALL THE SAMPLE NAMES FROM THE MERGED OBJECT

sample.names <- Final$File.Name
# Matrix data
# Loading feature table
feature_table <- read.delim("Metabolomics Datasets/M407_extracted_data_Joshua/c18neg/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable.txt") %>% select(-c("mz.min","mz.max",
                                                                                                                                                                      "NumPres.All.Samples", "NumPres.Biological.Samples",
                                                                                                                                                                      "median_CV", "PeakScore", "Qscore", "Max.Intensity",
                                                                                                                                                                      "NumPres.All.Samples" ,"NumPres.Biological.Samples"))
# Calculating min value except zero
# Value to exclude
# exclude_value <- 0
# 
# # Calculate the minimum value across all columns, excluding exclude_value
# min_value <- min(unlist(feature_table[feature_table != exclude_value]))
# 
# # Print the result
# print(min_value)

# Print the result
# print(min_except)
#Subsetting to only keep specific columns
feature_table$mztime <- paste0(feature_table$mz, "_", feature_table$time)

Final_df_new <- feature_table[, c("mztime", sample.names)]


# # Grouping and calculating median of redundant chem names
# new_Final_df <- Final_df_new %>%
#   group_by(Name) %>%
#   summarise(across(everything(), ~median(., na.rm = TRUE)))
# 

#converting df to matrix
mat <- as.matrix(Final_df_new)

#Converting first column to rownames
mat_new <- mat[,-1]
rownames(mat_new) <- mat[,1]


# check normalization
mat_numeric = apply(mat_new, 2, as.numeric)
rownames(mat_numeric)  = rownames(mat_new)

# boxplot(mat_numeric)

# Log transformation (added one cause of zero values)
# log_data <- log(mat_numeric+1)
log_data = log2(mat_numeric + 1)

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


# Converrting rownames to first column
scaled_df <- as.data.frame(scaled)
# adding mz column
scaled_df$mztime <- Final_df_new$mztime
rownames(scaled_df) <- NULL
# converting scaled to df

# Removing redundancy based on mz column by calculating median across every column
df <- scaled_df %>%
  group_by(mztime) %>%
  summarise(across(everything(), median))

df_t <- t(df)

#Making first row as colnames
colnames(df_t) <- df_t[1,]

# Deleting the first row which has now become the colname
df_t <- df_t[-1,]

# converting to a data frame 
data_f <- as.data.frame(df_t)
# Making rownames as a new column
data_f$sample <- rownames(data_f) 
# setting rownames as null
rownames(data_f) <- NULL

# setting sample as first column
# data_f <- as.numeric(data_f %>%  select(c("sample", everything())) )
data_f <- data_f %>%
  mutate_at(vars(-sample), as.numeric) %>% select(c("sample", everything()))


# extracting column names
metab_list <- colnames(data_f)[-1]

# Merging feature table and blood donor metadata and patient demographic

final_merged <- merge(data_f, Final, by.x = "sample", by.y = "File.Name")

# final_merged$blood_draw_to_irradiation <- as.numeric(final_merged$date_irradiation - final_merged$date_blood_draw)
# final_merged$irradiation_to_transfusion <- as.numeric(final_merged$trans_date - final_merged$date_irradiation)
# # final_merged$interaction_term <- as.numeric(final_merged$blood_draw_to_irradiation * final_merged$irradiation_to_transfusion)
# #checking the distribution of blood_draw_to_irradiation and irradiation to_transfusion
# final_merged <- mutate(final_merged,
#                        irradiation_to_transfusion_category = ifelse(irradiation_to_transfusion == 0, "zero", "morethanzero"))
# 
# # Now we factorize
# final_merged <- final_merged %>%
#   mutate(irradiation_to_transfusion_category = factor(irradiation_to_transfusion_category, levels = c("zero", "morethanzero")))

# As there is imbalance in the groups, mes response and no response, for the second model, we'll be taking oxygenation level after transfusion as well
# Also, we will add a column for age of blood and remove any blood more that 15 days old
# Additionally, factorize categorical variable

final_merged <- final_merged %>% mutate( age_of_blood= trans_date - date_blood_draw) %>% 
  filter(!age_of_blood > 15) %>% mutate(reci_sex= as.factor(sex.x)) %>% 
  mutate(mauc_mes_bef= log(mauc_mes_bef+1))  %>%
  mutate(mauc_mes_bef_scaled = scale(mauc_mes_bef)) %>% 
  mutate(oxy_change = mauc_mes_aft - mauc_mes_bef)

# checking the distribution of 
hist(final_merged$mauc_mes_bef)
hist(final_merged$mauc_mes_aft)
hist(final_merged$corrected_age)

# Now we perform glm
##
# Create an empty dataframe to store results
result_df <- data.frame(
  metabolite = character(),            # Column for metabolite name
  metab_beta = numeric(),     # Column for beta coefficient
  bo_beta = numeric(),        # Column for beta coefficient
  # age_of_blood = numeric(),
  # corrected_age = numeric(),
  p_value = numeric(),
  lower = numeric(),            # Column for lower confidence interval
  upper = numeric()            # Column for upper confidence interval
)

# Loop through metabolites
for (i in metab_list) {
  # Define the formula for the GLM
  formula <- as.formula(paste0("oxy_change   ~ `", i, "` + mauc_mes_bef"))
  
  # Fit the GLM
  model <- glm(formula, data = final_merged, family = gaussian()) 
  
  # Extract coefficients and confidence intervals
  temp <- summary(model)$coefficients
  conf_temp <- confint(model)
  
  # Extract p-value using ANOVA
  formula2 <- as.formula(paste0("oxy_change ~ 1"))
  model2 <- glm( formula2, data = final_merged)
  chunk <- anova(model, model2, test = "Chisq")
  p_value <- chunk$"Pr(>Chi)"[2]
  
  # Extract coefficients
  metab_beta <- temp[2, 1]
  bo_beta <- temp[3, 1]
  # age_of_blood <- temp[4, 1]
  # corrected_age <- temp[5, 1]
  
  # Extract confidence intervals
  lower_ <- conf_temp[2, 1]
  upper_ <- conf_temp[2, 2]
  
  # Add a row to result_df
  result_df <- rbind(result_df, c(i, metab_beta, bo_beta,  p_value, lower_, upper_))
}

# Assign column names to result_df
colnames(result_df) <- c("metabolite", "metab_beta", "bo_beta","p_value","lower", "upper")

# View the final dataframe
print(result_df)
str(result_df)
result_df$p_value <- as.numeric(result_df$p_value)
result_df$metab_beta <- as.numeric(result_df$metab_beta)
# Adjust p-values and calculate -log10
result_df$adjusted_p_value <- p.adjust(result_df$p_value, method = "BH")
result_df$log_adj_p_value <- -log10(as.numeric(result_df$adjusted_p_value))

## Checking normality of residual assumption
## loading necesssary 
library(tidyverse)
library(broom)
theme_set(theme_classic())

model.diag.metrics <- augment(model)

# Plotting
ggplot(model.diag.metrics, aes(mauc_mes_bef, mauc_mes_aft)) +
  geom_point() +
  stat_smooth(method = lm, se = FALSE) +
  geom_segment(aes(xend = mauc_mes_bef, yend = .fitted), color = "red", size = 0.3)

#
# Adjust margins and plot size
par(mar = c(3, 3, 2, 1))  # Adjust margins as needed
plot(model, mar = c(3, 3, 2, 1))  # Adjust plot margins

?par
########

result_df <- readRDS("glm_models/glm_base_model_c18.rds")

# Volcano plot
# add a column of NAs
result_df$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
result_df$diffexpressed[result_df$adjusted_p_value < 0.1 & result_df$metab_beta > 2] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
result_df$diffexpressed[result_df$adjusted_p_value< 0.1 & result_df$metab_beta < -2] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of mz features differentially expressed (NA in case they are not)
# result_df$delabel <- NA
# result_df$delabel[result_df$diffexpressed != "NO"] <- result_df$Name[result_df$diffexpressed != "NO"]

# de <- read.csv("test.csv")
# library(ggplot2)
# library(ggrepel)
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


ggplot(data=result_df, aes(x=metab_beta, y=-log10(p_value), col=diffexpressed)) +
  geom_point(size=0.5) + 
  theme_minimal() +
  scale_colour_manual(values = mycolors) +
  guides(color = guide_legend(title = NULL))

#saving the object
saveRDS(result_df, "glm_models/glm_oxy_change_base_model_c18.rds")

# Manhattan plot

# Add a new column to represent combinations of positive and negative values
result_df$combined_color <- case_when(
  result_df$metab_beta > 0  ~ "Both Positive",
  result_df$metab_beta < 0 ~ "Both Negative",
  TRUE ~ "Inconsistent Direction"
)

# Define colors for the different groups
color_palette_combined <- c("Both Positive" = "deeppink", "Both Negative" = "purple", "Inconsistent Direction" = "lightgreen")

significant_points <- result_df[result_df$log_p_value > -log10(0.01), ]

# Create the plot
manhattan_plot <- ggplot(result_df, aes(x = metab_beta, y = log_adj_p_value, color = combined_color)) +
  geom_point(size = 2, alpha= 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  labs(x = "", y = "-log10(adjusted p-value)", title = " ") +
  theme(axis.text.x = element_blank()) +  
  scale_color_manual(values = color_palette_combined, name = "Direction") +  # Change legend title for color
  scale_y_continuous(limits = c(0, 15)) +
  geom_text_repel(data = significant_points, aes(label = metabolite), vjust = -0.5, hjust = 1, size = 3, angle = 0, color= "black",
                  segment.color = "transparent")



# # CHECKKK
# i <- metab_list[1]f
# i is "100.0039_286.2"
# Define the formula for the GLM

# formula <- as.formula(paste0("`", i, "`", "~ blood_draw_to_irradiation + irradiation_to_transfusion + interaction_term"))
# Fit the GLM
# model <- glm(formula, data = final_merged, family = gaussian()) 

# model1 = metab ~ blood_draw_to_irradiation  + irradiation_to_tx + interaction_term
# model2 = metab ~ 1

#  model1 <- glm(`85.0295`~ blood_draw_to_irradiation + irradiation_to_transfusion + interaction_term, family = gaussian(), data = final_merged)
#  model2 <- glm(`85.0295` ~ 1, data = final_merged)
#  chunk <- anova(model1, model2, test = "Chisq")
#  temp <- summary(model1)$coefficients
# #
#  p_value <- chunk$"Pr(>Chi)"[2]
#  blood_draw_irr_beta <- temp[2, 1]
#  irr_trans_beta <- temp[3,1]
#  interaction_beta <- temp[4,1]
#  conf_temp <- confint(model1)
#  lower_ <- conf_temp[2, 1]
#  upper_ <- conf_temp[2, 2]

model1 <- glm(mauc_mes_aft ~ `177.0404_39.3` + reci_sex, family = gaussian(), data = final_merged)
sum <- summary(model1)
temp <- summary(model1)$coefficients
