library(data.table)
library(readxl)

## read UKB data
pas = readRDS("~/Downloads/UKB_other_vars_2025_v5.rds")

################################################
## overall code structure:
## extract the cancer types by icd10 code in the UK Biobank (incident and prevalent)
## extract the individuals with heterozygous or homozygous pLOF
## run logistic regression for each gene separately across all cancer types, adjusting for age/sex/10 genotype PCs
## make heatmap to show all results
## then make forest plot directly contrasting MSH2/3

## get the diseases:
library(dplyr)

# Extract all column names starting with f.41270 and f.41280
icd10_cols <- grep("^f\\.41270\\.0\\.", names(pas), value = TRUE)
date_cols <- grep("^f\\.41280\\.0\\.", names(pas), value = TRUE)

## get all ICD10
pi = pas[icd10_cols]
pi2 = as.data.frame(as.matrix(table(unlist(pi))))
pi2$ICD10 = rownames(pi2)

# Define a function to check DISEASE presence and find the earliest date
check_disease <- function(row) {
  # Extract ICD10 codes for the current row
  icd10_codes <- row[icd10_cols]
  
  # Check if any of the codes match the pattern for DISEASE (M45 or M45.x)
  # disease_indices <- grep("^M05", icd10_codes, value = FALSE) ## this is AS- M45 T1D- E10
  
  
  ## https://www.nature.com/articles/s41571-025-01054-6?fromPaywallRec=false Diaz Fig. 3- "spectrum"
  
  # disease_indices <- grep("^C18", icd10_codes, value = FALSE) ## CRC
  # disease_indices <- grep("^C54", icd10_codes, value = FALSE) ## endometrial (?)
  # disease_indices <- grep("^G10", icd10_codes, value = FALSE) ## Huntington's
  # disease_indices <- grep("^C56", icd10_codes, value = FALSE) ## ovarian
  # disease_indices <- grep("^C17", icd10_codes, value = FALSE) ## small bowel cancer
  # disease_indices <- grep("^C16", icd10_codes, value = FALSE) ## gastric cancer
  # disease_indices <- grep("^C61", icd10_codes, value = FALSE) ## prostate cancer
  # disease_indices <- grep("^C22", icd10_codes, value = FALSE) ## biliary cancer
  # disease_indices <- grep("^C71", icd10_codes, value = FALSE) ## GBM
  # disease_indices <- grep("^C67", icd10_codes, value = FALSE) ## bladder cancer
  # disease_indices <- grep("^C25", icd10_codes, value = FALSE) ## pancreatic cancer
  
  # disease_indices <- grep("^C44", icd10_codes, value = FALSE) ## cutaneous SCC
  # disease_indices <- grep("^C53", icd10_codes, value = FALSE) ## cervical
  # disease_indices <- grep("^C49", icd10_codes, value = FALSE) ## soft tissue
  # disease_indices <- grep("^C64", icd10_codes, value = FALSE) ## kidney
  # disease_indices <- grep("^C15", icd10_codes, value = FALSE) ## oesophageal
  # disease_indices <- grep("^C50", icd10_codes, value = FALSE) ## breast
  # disease_indices <- grep("^C43", icd10_codes, value = FALSE) ## head/neck (melanoma)
  # disease_indices <- grep("^C34", icd10_codes, value = FALSE) ## lung
  # disease_indices <- grep("^C73", icd10_codes, value = FALSE) ## thyroid
  
  # disease_indices <- grep("^C81", icd10_codes, value = FALSE) ## Hodgkin's lymphoma
  # disease_indices <- grep("^C859", icd10_codes, value = FALSE) ## Non-Hodgkin's lymphoma
  disease_indices <- grep("^C90", icd10_codes, value = FALSE) ## multiple myeloma
  
  if (length(disease_indices) > 0) {
    # Set has_disease to 1 if DISEASE is found
    has_disease <- 1
    
    # Extract corresponding dates for DISEASE occurrences
    disease_dates <- as.Date(row[date_cols[disease_indices]], format = "%Y-%m-%d")
    
    # Find the earliest non-NA date
    date_first <- min(disease_dates, na.rm = TRUE)
  } else {
    has_disease <- 0
    date_first <- NA
  }
  
  return(c(has_disease, date_first))
}

# Apply the function to each row of the dataframe
results <- t(apply(pas, 1, check_disease))

# Add the results to the dataframe
pas$has_disease <- results[, 1]
pas$date_first <- as.Date(results[, 2], origin = "1970-01-01")

## track if incident or prevalent case
pas$incident = ifelse(as.Date(pas$date_first) > as.Date(pas$f.53.0.0), 1, 0)

## incorporate assay date for antigen
pas$incident2 = ifelse(as.Date(pas$date_first) > as.Date(pas$f.53.0.0) & 
                         as.Date(pas$f.23048.0.0) < as.Date(pas$date_first), 1, 0)

table(pas$has_disease)

table(pas[pas$has_disease == 1,]$f.31.0.0)

# saveRDS(pas, "~/Downloads/pas_CRC.rds")
# saveRDS(pas, "~/Downloads/pas_endometrial_cancer.rds")
# saveRDS(pas, "~/Downloads/pas_ovarian.rds")
# saveRDS(pas, "~/Downloads/pas_smallbowel.rds")
# saveRDS(pas, "~/Downloads/pas_gastric.rds")
# saveRDS(pas, "~/Downloads/pas_prostate.rds")
# saveRDS(pas, "~/Downloads/pas_biliary.rds")
# saveRDS(pas, "~/Downloads/pas_GBM.rds")
# saveRDS(pas, "~/Downloads/pas_bladder.rds")
# saveRDS(pas, "~/Downloads/pas_pancreatic.rds")

# saveRDS(pas, "~/Downloads/pas_cutaneousSCC.rds")
# saveRDS(pas, "~/Downloads/pas_cervical.rds")
# saveRDS(pas, "~/Downloads/pas_softtissue.rds")
# saveRDS(pas, "~/Downloads/pas_kidney.rds")
# saveRDS(pas, "~/Downloads/pas_oesophageal.rds")
# saveRDS(pas, "~/Downloads/pas_breast.rds")
# saveRDS(pas, "~/Downloads/pas_headandneck.rds") ## this is actually melanoma
# saveRDS(pas, "~/Downloads/pas_lung.rds")
# saveRDS(pas, "~/Downloads/pas_thyroid.rds")

# saveRDS(pas, "~/Downloads/pas_hodgkin.rds")
# saveRDS(pas, "~/Downloads/pas_nonhodgkin.rds")
# saveRDS(pas, "~/Downloads/pas_multiplemyeloma.rds")


## load mutations

## maybe get muts in ALL genes

## or at least major drivers: BRCA1/2, TP53, APC, PTEN, CDH1, STK11, PALB2

## HC + LC (high confidence/low confidence)
# m = list.files("~/Downloads/MSH3_analysis", "raw", full.names = TRUE)

## HC only (high confidence)
m = list.files("~/Downloads/MSH3_analysis", "raw", full.names = TRUE)
m = m[grepl("HC_ONLY", m) == TRUE]

## function to obtain, for each gene, the individuals that have pLOF variants (het/hom/none) in that gene
mdf = data.frame()
for(j in m){
  
  print(j)
  
  gname = unlist(strsplit(basename(j), "_"))[1]
  
  print(gname)
  
  muts = as.data.frame(fread(j))
  
  for(i in colnames(muts)[7:ncol(muts)]){
    # print(i)
    # print(table(muts[,i]))
    print(which(colnames(muts) == i))
    tdf = data.frame(gname, i, paste(as.character(table(muts[,i])), collapse = ","), paste(muts[muts[,i] %in% c(0,1),]$FID, collapse = ","))
    colnames(tdf) = c("gene", "variant", "status", "IDs")
    mdf = rbind(mdf, tdf)
    # print('----------')
  }
}
mdf$nstat = as.numeric(apply(mdf, 1, function(x) unlist(strsplit(x["status"], ","))[1]))


# pas1 = readRDS("~/Downloads/pas_CRC.rds")
# pas2 = readRDS("~/Downloads/pas_endometrial_cancer.rds")

## all cancer types
ctypes = c("CRC", "endometrial_cancer", "ovarian", "smallbowel", 
           "gastric", "prostate", "biliary", "GBM", "bladder", 
           "pancreatic", "cutaneousSCC", "cervical", "softtissue", 
           "kidney", "oesophageal", "breast", "headandneck", "lung", 
           "thyroid", "hodgkin", "nonhodgkin", "multiplemyeloma")

## cancer types in men (coded S1) only
ctypes = c("CRC", "smallbowel", 
           "gastric", "prostate", "biliary", "GBM", "bladder", 
           "pancreatic", "cutaneousSCC", "softtissue", 
           "kidney", "oesophageal", "headandneck", "lung", 
           "thyroid", "hodgkin", "nonhodgkin", "multiplemyeloma")

## cancer types in women (coded S0) only
ctypes = c("CRC", "endometrial_cancer", "ovarian", "smallbowel", 
           "gastric", "biliary", "GBM", "bladder", 
           "pancreatic", "cutaneousSCC", "cervical", "softtissue", 
           "kidney", "oesophageal", "breast", "headandneck", "lung", 
           "thyroid", "hodgkin", "nonhodgkin", "multiplemyeloma")

## logistic regression analysis
## the key variable is "mut"- which is 1 if the individual has a het OR hom pLOF variant (to boost power, I grouped them together)
## 0 otherwise
## then adjust for age, sex, 10 genotype PCs
## debated including smoking status, but smoking does not predict cancer equally across all types, so I left out for now.
mmf = data.frame()
for(i in ctypes){
  print(i)
  pas = readRDS(sprintf("~/Downloads/pas_%s.rds", i))
  pas$sex = paste("S", pas$f.31.0.0, sep = "")
  
  for(j in unique(mdf$gene)){
    print(j)
    pas$mut = ifelse(pas$f.eid %in% unique(unlist(strsplit(mdf[mdf$gene == j,]$IDs, ","))), 1, 0)
    
    ## everyone
    model <- glm(has_disease ~ mut + f.21022.0.0 + sex + f.22009.0.1 + f.22009.0.2 +
                   f.22009.0.3 + f.22009.0.4 + f.22009.0.5 +
                   f.22009.0.6 + f.22009.0.7 + f.22009.0.8 +
                   f.22009.0.9 + f.22009.0.10, data = pas, family = binomial)
    
    ## women only
    # model <- glm(has_disease ~ mut + f.21022.0.0 + f.22009.0.1 + f.22009.0.2 +
    #                f.22009.0.3 + f.22009.0.4 + f.22009.0.5 +
    #                f.22009.0.6 + f.22009.0.7 + f.22009.0.8 +
    #                f.22009.0.9 + f.22009.0.10, data = pas[pas$sex == "S0",], family = binomial)
    
    ## men only
    # model <- glm(has_disease ~ mut + f.21022.0.0 + f.22009.0.1 + f.22009.0.2 +
    #                f.22009.0.3 + f.22009.0.4 + f.22009.0.5 +
    #                f.22009.0.6 + f.22009.0.7 + f.22009.0.8 +
    #                f.22009.0.9 + f.22009.0.10, data = pas[pas$sex == "S1",], family = binomial)
    
    # Summary of the model
    summary_model <- summary(model)
    
    # Extract coefficients
    coefficients <- summary_model$coefficients
    
    # Extract beta for the specified allele
    beta <- coefficients["mut", "Estimate"]
    
    # Calculate Odds Ratio
    odds_ratio <- exp(beta)
    
    beta_CI_low = confint.default(model)[2,][1]
    beta_CI_high = confint.default(model)[2,][2]
    
    OR_CI_low = exp(beta_CI_low)
    OR_CI_high = exp(beta_CI_high)
    
    pval = coefficients[2,4]
    
    ## everyone
    tdf = data.frame(i, table(pas$has_disease)[2], table(pas$has_disease)[1], j, table(pas$mut)[2], table(pas$mut)[1],
                     pval, beta, beta_CI_low, beta_CI_high, odds_ratio, OR_CI_low, OR_CI_high)
    
    ## women only
    # tdf = data.frame(i, table(pas[pas$sex == "S0",]$has_disease)[2], 
    #                  table(pas[pas$sex == "S0",]$has_disease)[1], j, 
    #                  table(pas[pas$sex == "S0",]$mut)[2], 
    #                  table(pas[pas$sex == "S0",]$mut)[1],
    #                  pval, beta, beta_CI_low, beta_CI_high, odds_ratio, OR_CI_low, OR_CI_high)
    
    ## men only
    # tdf = data.frame(i, table(pas[pas$sex == "S1",]$has_disease)[2],
    #                  table(pas[pas$sex == "S1",]$has_disease)[1], j,
    #                  table(pas[pas$sex == "S1",]$mut)[2],
    #                  table(pas[pas$sex == "S1",]$mut)[1],
    #                  pval, beta, beta_CI_low, beta_CI_high, odds_ratio, OR_CI_low, OR_CI_high)
    
    colnames(tdf) = c("cancer", "Ncase", "Ncontrol", "gene", "Nmut", "Nnomut", "pval", "beta", "beta_CI_low", "beta_CI_high", "OR", "OR_CI_low", "OR_CI_high")
    mmf = rbind(mmf, tdf)
  }
  print('-------------------')
}  

# saveRDS(mmf, "~/Downloads/mmf_HC_ONLY.rds")
# saveRDS(mmf, "~/Downloads/mmf_HC_LC.rds")

# saveRDS(mmf, "~/Downloads/mmf_HC_ONLY_women.rds")
# saveRDS(mmf, "~/Downloads/mmf_HC_ONLY_men.rds")

# mmf = readRDS("~/Downloads/mmf_HC_LC.rds")
# mmf = readRDS("~/Downloads/mmf_HC_ONLY.rds")
mmf = readRDS("~/Downloads/mmf_HC_ONLY.rds")


f4 = matrix(0, nrow = length(unique(mmf$gene)), ncol = length(unique(mmf$cancer)))
rownames(f4) = sort(unique(mmf$gene))
colnames(f4) = sort(unique(mmf$cancer))

for(i in rownames(f4)){
  print(i)
  temp = mmf[mmf$gene == i,]
  for(j in colnames(f4)){
    temp2 = temp[temp$cancer == j,]
    
    ## Bonferroni
    if(temp2$pval <= 0.1/nrow(mmf)){f4[i,j] = temp2$OR}
    if(temp2$pval > 0.1/nrow(mmf)){f4[i,j] = NA}
    
    ## Unadjusted
    # if(temp2$pval <= 0.05){f4[i,j] = temp2$OR}
    # if(temp2$pval > 0.05){f4[i,j] = NA}
  }
}
f4[is.na(f4)] = 0

f4 = f4[,c("biliary", "bladder", "CRC", "endometrial_cancer", "gastric", 
           "GBM", "ovarian", "pancreatic", "prostate", "smallbowel", ## Lynch solid
           "breast", "cervical", "cutaneousSCC", "kidney", "lung", 
           "headandneck", "oesophageal", "softtissue", "thyroid", ## Non-Lynch solid
           "hodgkin", "nonhodgkin", "multiplemyeloma")] ## blood

# Create a custom color palette
# Start with white for 0, then use shades of red for higher values
reds <- brewer.pal(9, "Reds")
custom_colors <- c("white", reds)

# Define breaks so that 0 maps to white
# You need to set breaks manually to ensure 0 is treated specially
breaks <- seq(min(f4), max(f4), length.out = length(custom_colors))

display_numbers <- ifelse(f4 != 0, round(f4,2), "")

for(i in rownames(f4)){
  rownames(f4)[rownames(f4) == i] = paste(c(i, " ", "(", "N = ", unique(mmf[mmf$gene == i,]$Nmut), ")"), collapse = "")
}

for(i in colnames(f4)){
  colnames(f4)[colnames(f4) == i] = paste(c(i, " ", "(", "N = ", unique(mmf[mmf$cancer == i,]$Ncase), ")"), collapse = "")
}

colnames(f4)[grepl("CRC", colnames(f4)) == TRUE] = gsub("CRC", "colorectal", colnames(f4)[grepl("CRC", colnames(f4)) == TRUE])
colnames(f4)[grepl("cutaneousSCC", colnames(f4)) == TRUE] = gsub("cutaneousSCC", "cutaneous_SCC", colnames(f4)[grepl("cutaneousSCC", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("endometrial_cancer", colnames(f4)) == TRUE] = gsub("endometrial_cancer", "endometrial", colnames(f4)[grepl("endometrial_cancer", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("GBM", colnames(f4)) == TRUE] = gsub("GBM", "glioblastoma", colnames(f4)[grepl("GBM", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("headandneck", colnames(f4)) == TRUE] = gsub("headandneck", "melanoma", colnames(f4)[grepl("headandneck", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("hodgkin", colnames(f4)) == TRUE & grepl("non", colnames(f4)) == FALSE] = gsub("hodgkin", "hodgkins_lymphoma", colnames(f4)[grepl("hodgkin", colnames(f4)) == TRUE & grepl("non", colnames(f4)) == FALSE]) 
colnames(f4)[grepl("multiplemyeloma", colnames(f4)) == TRUE] = gsub("multiplemyeloma", "multiple_myeloma", colnames(f4)[grepl("multiplemyeloma", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("nonhodgkin", colnames(f4)) == TRUE] = gsub("nonhodgkin", "nonhodgkins_lymphoma", colnames(f4)[grepl("nonhodgkin", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("smallbowel", colnames(f4)) == TRUE] = gsub("smallbowel", "small_bowel", colnames(f4)[grepl("smallbowel", colnames(f4)) == TRUE]) 
colnames(f4)[grepl("softtissue", colnames(f4)) == TRUE] = gsub("softtissue", "soft_tissue", colnames(f4)[grepl("softtissue", colnames(f4)) == TRUE]) 


pheatmap(f4, show_rownames = TRUE, show_colnames = TRUE, scale = "none",
         cluster_rows = FALSE, cluster_cols = FALSE, 
         border_color = "black", na_col = "white",
         display_numbers = display_numbers,
         number_color = "black",
         color = custom_colors, breaks = breaks, 
         angle_col = 315, fontsize_col = 9, fontsize_row = 9, fontsize_number = 8) #brewer.pal(10, name = "Reds")


## forest plot for MSH3 vs MSH2
## maybe put them both on the same plot
m2 = mmf[mmf$gene %in% c("MSH3", "MSH2") & (! mmf$cancer %in% c("smallbowel", "hodgkin", "thyroid", 
                                                                "softtissue", "GBM", "nonhodgkin", 
                                                                "pancreatic", "multiplemyeloma", "ovarian", 
                                                                "headandneck", "cervical")),]
m2$color = ifelse(m2$pval <= 0.05, "red", "black")
m2$gene = factor(m2$gene, levels = rev(c("MSH2", "MSH3")))

m3 = m2[m2$gene == "MSH2",]
m2$cancer = factor(m2$cancer, levels = m3[order(m3$OR),]$cancer)

ggplot(data=m2, aes(y=cancer, x=OR, xmin=OR_CI_low, xmax=OR_CI_high, group = gene, shape = gene)) +
  # geom_pointrange(size = 1, color = m2$color) + 
  
  geom_pointrange(color = m2$color, position = position_dodge(width = 1), size = 0.8) + 
  
  geom_vline(xintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  # coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + 
  ylab("") +
  theme_bw() +
  # ylim(c(-9e-05,0.0002)) +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 0, hjust=0.7)) +
  theme(axis.text.y = element_text(size = 14, color = "black")) +
  theme(axis.title.x = element_text(size = 14, color = "black")) +
  theme(axis.title.y = element_text(size = 14, color = "black")) +
  scale_shape_manual(values = c(16, 17)) +
  guides(shape = guide_legend(reverse = TRUE))
