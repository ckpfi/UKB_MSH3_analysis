library(data.table)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

x = fread(sprintf("UKB_450k_exome_Broad_chr%s_annot.txt.gz", args[2]))
x2 = x[x$Gene_symbol == args[1],] ## MSH3 (5) ## MLH1 (3) ## MSH2 (2) ## MSH6 (2) ## PMS2 (7)
f = fread("/hpfs/projects/ukbiobank/genotypes/UKB_450k_exome_Broad/UKB_450k_exome_Broad_EUR_QC_freq.afreq")
x3 = x2[x2$Variant %in% f$ID,]
table(x3$LoF, x3$ENST)
length(unique(x3[x3$LoF == "HC",]$Variant))
#length(unique(x3[x3$LoF %in% c("HC", "LC"),]$Variant))
#write.table(unique(x3[x3$LoF %in% c("HC", "LC"),]$Variant), "PMS1_HC_LC_LOF.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

write.table(unique(x3[x3$LoF %in% c("HC"),]$Variant), sprintf("%s_HC_ONLY_LOF.txt", args[1]), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
