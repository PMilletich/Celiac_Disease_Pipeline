library(ggplot2)
library(ggforce)
library(vegan)
library(ggpubr)
library(DESeq2)
library(phyloseq)
library(dplyr)
'%notin%' = Negate('%in%')

###### 
#File import 
setwd("~/Desktop/Reviewer_Edits_V2/")
taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
Total_reads = readRDS("Filters_PS.RDS")
ps.RA.ALL = transform_sample_counts(Total_reads, function(x) x / sum(x) )
OTU_All = data.frame(otu_table(ps.RA.ALL))

list = c()
for (i in 1:ncol(OTU_All)) {
  current_list = OTU_All[,i]
  list = c(list, current_list)
}

list= subset(list, list != 0)
summary(list)
 
Sample_data = read.csv("Celiac_Samples.csv")

Rel_int = 1e6
alpha = 0.05
Control_Num = 2
seed_number = 1
seed_max = 100

current_glom_1 = "Genus"

current_method = "Matched"
Actual_sign = data.frame()
for (seed_number in 1:seed_max) {
  if (seed_number%%10 == 0 ) {
    print(seed_number)
  }
  ps.filtered = Total_reads
  #####################
  #Sample Matching with n Controls
  subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
  subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
  
  final_subset = subset_cases
  Subset_1 = c("Region","Siblings_at_birth")
  
  Matched = list()
  Count_1 = 0
  current_id = "19289"#unique(subset_df_1$ID)[1]
  for (current_id in unique(subset_cases$ID)){
    current_row = subset(subset_cases, subset_cases$ID == current_id)
    current_row_subset = current_row[,Subset_1]
    current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
    
    if (current_id == "2061") {
      current_row_subset = data.frame("Region"= current_row_subset)
    }
    
    Matched_subset = merge(subset_controls, current_row_subset)
    Matched = c(Matched, nrow(Matched_subset))
    # print(paste(current_id, nrow(Matched_subset)))
    
    if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
      set.seed(seed_number)
      Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
      Count_1 = Count_1 + 1
    } else {
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
    }
    subset_controls = subset(subset_controls, ! subset_controls$ID %in% final_subset$ID)
  }
  
  
  final_subset = unique(final_subset)
  table(final_subset$Autoimmune_Disorder)
  
  rownames(final_subset) = final_subset$ID
  sample_data(ps.filtered) = sample_data(final_subset)
  
  #################################################################
  #Relative Abundance; https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/transform_sample_counts
  ps.RA = transform_sample_counts(ps.filtered, function(x) x / sum(x) )
  
  #Taxa Identification 
  row.names(taxa_identifiers) = taxa_identifiers$ASV
  taxa_identifiers[,current_glom_1] = gsub("_V2", "", taxa_identifiers[,current_glom_1])
  taxa_identifiers[,current_glom_1] = gsub("_ASV", ".", taxa_identifiers[,current_glom_1])
  taxa_identifiers = taxa_identifiers[,c("ASV", current_glom_1)]
  
  ps.current = ps.RA
  
  ###########
  #Data Fluffing
  count_table = data.frame(t(otu_table(ps.current)))
  count_table= round(count_table * Rel_int)
  count_table$ASV = row.names(count_table)
  count_table = merge(count_table, taxa_identifiers)
  row.names(count_table) = count_table[,current_glom_1]
  Identifier_list = unique(rownames(count_table))
  tax_table_all_group = unique(count_table[,current_glom_1])
  count_table$ASV = NULL
  count_table[,current_glom_1]= NULL
  count_table = data.frame(t(count_table))
  row.names(count_table) = gsub("X", "", row.names(count_table))
  count_table$ID = row.names(count_table)
  count_table_2 = merge(count_table, final_subset)
  AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
  Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
  count_table_2$Autoimmune_2_groups = factor(count_table_2$Autoimmune_2_groups, 
                                             levels = unique(count_table_2$Autoimmune_2_groups))
  
  OTU = data.frame(otu_table(ps.current))
  #print(mean(as.matrix(OTU)))
  OTU =  round(OTU * Rel_int)
  OTU[is.na(OTU)] = 0
  OTU <- mutate_all(OTU, function(x) as.integer(as.character(x)))
  otu_table(ps.current) = otu_table(OTU, taxa_are_rows = F)
  diagdds = suppressWarnings(phyloseq_to_deseq2(ps.current, ~ Autoimmune_2_groups))
  diagdds = estimateSizeFactors(diagdds, type = "poscounts")
  
  diagdds = DESeq(diagdds, test="Wald", fitType="local", quiet = TRUE)
  resultsNames(diagdds)
  res = results(diagdds, cooksCutoff = FALSE)
  
  sigtab = data.frame(res[which(res$padj < alpha), ])
  
  if (nrow(sigtab)> 0) {
    sigtab = res[which(res$padj < alpha), ]
    sigtab = data.frame(sigtab)
    sigtab$ASV = rownames(sigtab)
    sigtab <- sigtab[,colSums(is.na(sigtab))<nrow(sigtab)]
    sigtab = merge(sigtab, taxa_identifiers)
    sigtab[,current_glom_1] = gsub("-", ".", sigtab[,current_glom_1])
    sigtab[,current_glom_1] = gsub("/", ".", sigtab[,current_glom_1])
    
    Identifier_list = unique(sigtab[,current_glom_1])
    current_taxa = Identifier_list[1]
    significant_data = data.frame()
    for (current_taxa in Identifier_list) {
      AI[,current_taxa] = as.numeric(AI[,current_taxa])
      Control[,current_taxa] = as.numeric(Control[,current_taxa])
      AI_prev = round(nrow(subset(AI, AI[,current_taxa] > 0))/nrow(AI), 3)*100
      Control_prev = round(nrow(subset(Control, Control[,current_taxa] > 0))/nrow(Control), 3)*100
      count_table_2[,current_taxa] = as.numeric(count_table_2[,current_taxa])
      count_table_2$Genera  = count_table_2[,current_taxa]
      P_value = subset(sigtab, sigtab[,current_glom_1] == current_taxa)
      P_value = P_value$padj
      current_row = data.frame("Taxa" = current_taxa, "AI_prev" = AI_prev, 
                               "Control_prev" = Control_prev, "P_val" = P_value)
      significant_data = rbind(significant_data,current_row )
    }
    sigtab = sigtab[,c(current_glom_1, "log2FoldChange", "baseMean")]
    colnames(sigtab) = c("Taxa_Name", "log2FoldChange", "baseMean" )
    Actually_significant = significant_data
    Actually_significant$Taxa_Name = Actually_significant$Taxa
    Actually_significant$Taxa = current_glom_1
    Actually_significant_1 = merge(Actually_significant, sigtab)
    Actual_sign = rbind(Actual_sign, Actually_significant_1)
  }
}

print(table(Actual_sign$Taxa_Name))
write.csv(Actual_sign, paste("./CSV_Output/Deseq_", current_method, "RelAbun_2.csv", sep = ""), row.names = F)
