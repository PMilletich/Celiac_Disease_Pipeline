library(DESeq2)
library(dplyr) #mutate_all
'%notin%' = Negate('%in%')

################################################################# 
#File import 
#################################################################
#File with 16s Sequence and taxonomical classification 
taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
row.names(taxa_identifiers) = taxa_identifiers$ASV
taxa_identifiers = taxa_identifiers[,c("ASV", "Genus")]

#From Step 1 
Total_reads = readRDS("Filters_PS.RDS")
Sample_data = read.csv("Celiac_Samples.csv")

#################################################################
#Default Values 
#################################################################
alpha = 0.05 #Significance threshold 
Control_Num = 2 #Numnber of Cases per control 
seed_max = 100 #Max number of iterations 

#Create dataframe to append to 
Actual_sign = data.frame()

for (seed_number in 1:seed_max) {
  #Print 10s to track progress 
  if (seed_number%%10 == 0 ) {
    print(seed_number)
  }
  #Create current iterations phyloseq object 
  ps.filtered = Total_reads
  #####################
  #Sample Matching with n Controls
  subset_cases = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
  subset_controls = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")
  
  #Save Cases to final dataframe
  final_subset = subset_cases
  #Create list of matching variables: Step 2
  Subset_1 = c("Region","Siblings_at_birth")
  
  #Create list of number of controls for analysis 
  Matched = list()
  Count_1 = 0
  
  #For all Celiac cases IDs
  for (current_id in unique(subset_cases$ID)){
    #Subset dataframe for ID, and only metadata columns of interest 
    current_row = subset(subset_cases, subset_cases$ID == current_id)
    current_row_subset = current_row[,Subset_1]
    #Remove columns with NA
    current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
    
    #Case ID 2061 only had Region listed, create a dataframe accordingly 
    if (current_id == "2061") {
      current_row_subset = data.frame("Region"= current_row_subset)
    }
    
    #Merge controls matching columns of interest 
    Matched_subset = merge(subset_controls, current_row_subset)
    #Append number of available controls to Matched list for further analysis 
    Matched = c(Matched, nrow(Matched_subset))
    
    #Subset the controls based on max number of controls from line ~18
    if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
      set.seed(seed_number) #Set seed based on iteration number 
      Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
      Count_1 = Count_1 + 1
    } else {
      set.seed(seed_number)
      Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
      final_subset = rbind(final_subset, Matched_subset)
    }
    #Remove subset controls from the available control list 
    subset_controls = subset(subset_controls, ! subset_controls$ID %in% final_subset$ID)
  }
  
  #Add select case-control cohort to the phyloseq object 
  rownames(final_subset) = final_subset$ID
  sample_data(ps.filtered) = sample_data(final_subset)
  
  #################################################################
  #Relative Abundance; https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/transform_sample_counts
  ps.RA = transform_sample_counts(ps.filtered, function(x) x / sum(x) )
  
  #################################################################
  #qPCR; 16s rRNA reads/g
  qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
  row.names(qpcr_data) = qpcr_data$ID
  qpcr_data$ID = NULL
  
  #Merge the Relative Abundance table and the qpcr table 
  otu_RA = data.frame(otu_table(ps.RA))
  otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
  row.names(otu_RA_1) = otu_RA_1$Row.names
  otu_RA_1$Row.names = NULL

  #Multiply all the Relabundance columns by the reads/g column
  for(i in 2:length(colnames(otu_RA_1))) {
    otu_RA_1[is.na(otu_RA_1[,i]),i] = 0
    otu_RA_1[,i] <- suppressWarnings(as.integer(round(otu_RA_1[,1] * otu_RA_1[, i])))
  }
  ps.RA.qpcr_1 = ps.RA
  otu_table(ps.RA.qpcr_1) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
  
  #use the total abundance (qpcr) data for all following analysis
  ps.current = ps.RA.qpcr_1
                                  
  ###########
  #Create Dataframe for analysis
  #Create data frame of OTU Table and merge with Taxa_identifiers to find Genus
  count_table = data.frame(t(otu_table(ps.current)))
  count_table$ASV = row.names(count_table)
  count_table = merge(count_table, taxa_identifiers)
  row.names(count_table) = count_table[,current_glom_1]
  #Remove old taxa_identifier columns                                 
  count_table$ASV = NULL
  count_table[,current_glom_1]= NULL
  #Transpose (IDs as rows), and merge with metadata                                
  count_table = data.frame(t(count_table))
  row.names(count_table) = gsub("X", "", row.names(count_table))
  count_table$ID = row.names(count_table)
  count_table_2 = merge(count_table, final_subset)
  #Create subset of Celiac and controls                                
  AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
  Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
  count_table_2$Autoimmune_2_groups = factor(count_table_2$Autoimmune_2_groups, 
                                             levels = unique(count_table_2$Autoimmune_2_groups))
  #Replace any NA's and make sure values are integers 
  OTU = data.frame(otu_table(ps.current))
  OTU[is.na(OTU)] = 0
  OTU <- mutate_all(OTU, function(x) as.integer(as.character(x)))
  otu_table(ps.current) = otu_table(OTU, taxa_are_rows = F)
  
  #Run DESEQ; https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/estimateSizeFactors
  diagdds = suppressWarnings(phyloseq_to_deseq2(ps.current, ~ Autoimmune_2_groups))
  diagdds = estimateSizeFactors(diagdds, type = "poscounts")
  diagdds = DESeq(diagdds, test="Wald", fitType="local", quiet = TRUE)
  resultsNames(diagdds)
  res = results(diagdds, cooksCutoff = FALSE)
  #Subset based on significance threshold; line ~20
  sigtab = data.frame(res[which(res$padj < alpha), ])
  
  #If there are significant taxa 
  if (nrow(sigtab)> 0) {
    sigtab$ASV = rownames(sigtab)
    #Remove empty columns 
    sigtab <- sigtab[,colSums(is.na(sigtab))<nrow(sigtab)]
    #Merge with taxa_identifiers 
    sigtab = merge(sigtab, taxa_identifiers)
    #Replace special characters 
    sigtab[,current_glom_1] = gsub("-", ".", sigtab[,current_glom_1])
    sigtab[,current_glom_1] = gsub("/", ".", sigtab[,current_glom_1])
    
    #Create list of significant taxa 
    Identifier_list = unique(sigtab$Genus)
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
      #Save the taxa name, prevalence for each group, and padj value into dataframe
      current_row = data.frame("Taxa" = current_taxa, "AI_prev" = AI_prev, 
                               "Control_prev" = Control_prev, "P_val" = P_value)
      significant_data = rbind(significant_data,current_row )
    }
  }
  sigtab = sigtab[,c(current_glom_1, "log2FoldChange", "baseMean")]
  colnames(sigtab) = c("Taxa_Name", "log2FoldChange", "baseMean" )
  Actually_significant = significant_data
  Actually_significant$Taxa_Name = Actually_significant$Taxa
  Actually_significant$Taxa = current_glom_1
  Actually_significant_1 = merge(Actually_significant, sigtab)
  Actual_sign = rbind(Actual_sign, Actually_significant_1)
}

#Save dataframe as CSV
write.csv(Actual_sign, paste("./CSV_Output/Deseq_", current_method, "_2.csv", sep = ""), row.names = F)
