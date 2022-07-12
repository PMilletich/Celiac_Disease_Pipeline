#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr)

#################################################################
#File Import 
physeq_f = readRDS("Filters_PS.RDS")
Sample_data = data.frame(sample_data(physeq_f))

taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
row.names(taxa_identifiers) = taxa_identifiers$ASV
taxa_identifiers = taxa_identifiers[,c("ASV", "Genus")]


#################################################################
#Default Values 
setwd("~/Desktop/ABIS/Celiac_Paper/Frontiers_Jun13/Reviewer_Edits_v2/")
Control_Num = 2
seed_max = 100 
seed_number = 1 
current_method = "Reads"
iteration_threshold = 50 

#Create list of matching variables: Step 2
Subset_1 = c("Region","Siblings_at_birth")

subset_cases_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

#################################################################
#Relative Abundance 
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

#################################################################
#qPCR; Total Abundance 
qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
rownames(qpcr_data) = qpcr_data$ID
qpcr_data$ID = NULL

#Merge the Relative Abundance table and the qpcr table 
otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL

for(i in 2:length(names(otu_RA_1))) {
  otu_RA_1[,i] <- round(otu_RA_1[,1] * otu_RA_1[, i])
}

otu_RA_1$copies_16s_per_gram_stool = NULL

ps.RA.qpcr = ps.RA
otu_table(ps.RA.qpcr) = otu_table(otu_RA_1, taxa_are_rows = FALSE)
#tax_table(ps.RA) = tax_table(ps.RA.qpcr)

for (current_method in c("Reads", "RelativeAbundance")) {
  total_data = data.frame()
  best_df = data.frame()
  for (seed_number in 1:seed_max) {
    iteration = seed_number
    
    #Sample Matching with n Controls
    subset_cases = subset_cases_OG
    subset_controls = subset_controls_OG
    
    #Save Cases to final dataframe
    PIME_data = subset_cases
    PIME_data$Group = "fCD"
    
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
      Matched_subset$Group = "Controls"
      
      #Subset the controls based on max number of controls from line ~20
      if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
        set.seed(seed_number) #Set seed based on iteration number 
        Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
        PIME_data = rbind(PIME_data, Matched_subset)
      } else { #If there are no matched controls, randomly select you
        set.seed(seed_number)
        Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
        PIME_data = rbind(PIME_data, Matched_subset)
      }
      #Remove subset controls from the available control list 
      subset_controls = subset(subset_controls, ! subset_controls$ID %in% PIME_data$ID)
      
    }
    table(PIME_data$Group)
    
    #Load Relative or Total Abundance phyloseq object 
    if (current_method == "Reads") {
      physeq_f = ps.RA.qpcr
    } else {
      physeq_f = ps.RA
    }
    #Add select case-control cohort to the phyloseq object 
    rownames(PIME_data) = PIME_data$ID
    sample_data(physeq_f) = sample_data(PIME_data)
    physeq_f #[ 175 taxa and 78 samples ]
    
    ################################################
    #PIME
    print(pime.oob.error(physeq_f, "Group"))
    per_variable_obj= pime.split.by.variable(physeq_f, "Group")
    
    prevalences=pime.prevalence(per_variable_obj)
    set.seed(42)

    best.prev=pime.best.prevalence(prevalences, "Group")

    #Add information for the first 10 iterations 
    if (seed_number <= 10) {
      print(paste(seed_number, sum(phyloseq::sample_sums(physeq_f))))
      
      best.df = data.frame(best.prev$`OOB error`)
      best.df$Iteration = seed_number
      if (seed_number == 1) {
        best_df= best.df
      } else {
        best_df = cbind(best_df, best.df)
      }
    }
  
    imp60=best.prev$`Importance`$`Prevalence 60`
    prevalence.60 = prevalences$`60`
    input_ord = ordinate(prevalence.60, "PCoA" , "binomial")
    
    if (iteration%%25 == 0) {
      #Create PCOA graph
      graph = plot_ordination(prevalence.60, input_ord , color = "Group")+
        stat_ellipse(aes(group = Group)) + 
        ggtitle(paste("Iteration: ", iteration, sep = ""))
      assign(paste("graph_", iteration, sep = ""), graph)
    }
    
    imp60_subset = subset(imp60, imp60$MeanDecreaseAccuracy > 0 & 
                            imp60$Controls > 0 & 
                            imp60$fCD > 0)
    
    ###########
    #Create Dataframe for analysis
    #Create data frame of OTU Table and merge with Taxa_identifiers to find Genus
    count_table = data.frame(t(otu_table(prevalence.60)))
    count_table$ASV = row.names(count_table)
    count_table = merge(count_table, taxa_identifiers)
    row.names(count_table) = count_table$Genus
    #Remove old taxa_identifier columns                                 
    count_table$ASV = NULL
    count_table$Genus= NULL
    #Transpose (IDs as rows), and merge with metadata 
    count_table = data.frame(t(count_table))
    row.names(count_table) = gsub("X", "", row.names(count_table))
    count_table$ID = row.names(count_table)
    count_table_2 = merge(count_table, PIME_data)
    #Create subset of Celiac and controls                                
    AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
    Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
    
    #Find Mean Abundance for all selected Taxa 
    mean_data = data.frame()
    for (current_genus in unique(imp60_subset$Genus)) {
      current_subset = subset(imp60_subset, imp60_subset$Genus == current_genus)
      current_subset = current_subset[,c("Genus", "Controls", "fCD", "MeanDecreaseAccuracy")]
      current_genus = gsub("/", ".", current_genus)
      AI_mean = mean(as.numeric(AI[,current_genus]))
      Control_mean = mean(as.numeric(Control[,current_genus]))
      current_subset$AI_Mean = AI_mean
      current_subset$Control_Mean = Control_mean
      mean_data= rbind(mean_data, current_subset)
    }
    #Bind all the averaged data 
    total_data = rbind(total_data, mean_data)
  }
  #Graph sample iterations 
  final_graph = ggarrange(graph_25, graph_50, 
                          graph_75, graph_100, 
                          ncol = 2, nrow = 2, common.legend = T)
  final_graph = annotate_figure(final_graph, current_method)
  
  #Assign and save final dataframe and graph
  assign(paste(current_method, "total", sep = ""), total_data)
  assign(paste("PCOA", current_method, sep = "_"), final_graph)
}
#Merge together both PCOA graphs
final_graph = ggarrange(`PCOA_Reads`, PCOA_RelativeAbundance)

#Save final merged graph 
jpeg("PIME_60.jpeg", res = 400, height = 2000, width = 4000)
final_graph
dev.off()

#Create bar graphs of averages 
for (current_method in c("Reads/g", "Relative Abundance")) {
  #Load Averages of each iteration based on method
  if (current_method == "Reads/g") {
    total_data = Readstotal
  } else {
    total_data = RelativeAbundancetotal 
  }
  
  #Count number of counts by iterations and remove by iteration threshold ~line 24
  count_data = data.frame(table(total_data$Genus))
  count_data = subset(count_data, count_data$Freq >= iteration_threshold)
  #Only keep those that are above the threshold
  total_data_subset = subset(total_data, 
                             total_data$Genus %in% count_data$Var1)

  #Create dataframes of just fCD or Controls 
  total_data_fCD = total_data_subset[,c("AI_Mean", "MeanDecreaseAccuracy", "Genus")]
  colnames (total_data_fCD) = c("Value", "Accuracy", "Genus")
  total_data_fCD$Group = "fCD"
  
  total_data_Control = total_data_subset[,c("Control_Mean", "MeanDecreaseAccuracy", "Genus")]
  colnames (total_data_Control) = c("Value", "Accuracy", "Genus")
  total_data_Control$Group = "Control"
  
  fCD_data = data.frame()
  Control_data = data.frame()
  #For all the subsetted Genera...
  for (current_genus in unique(total_data_subset$Genus)) {
    AI_subset = subset(total_data_fCD, total_data_fCD$Genus == current_genus)
    AI_mean = mean(AI_subset$Value)
    
    Control_subset = subset(total_data_Control, total_data_Control$Genus == current_genus)
    Control_mean = mean(Control_subset$Value)
    #Separate by Trend and create graphs
    if (AI_mean > Control_mean) {
      fCD_data = rbind(fCD_data,AI_subset)
      fCD_data = rbind(fCD_data,Control_subset)
    } else if (AI_mean == Control_mean) {
      print(paste("Uh oh", current_genus))
      } else {
      Control_data = rbind(Control_data,AI_subset)
      Control_data = rbind(Control_data,Control_subset)
    }
  }
  
  fCD_graph = ggplot(fCD_data,  aes(x = Genus, y = Value, color = Group)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank()) + ylab(current_method) +
    ggtitle(paste(current_method, "\nGreater in fCD")) 
  
  Control_graph = ggplot(Control_data,  aes(x = Genus, y = Value, color = Group)) + 
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank()) + ylab(current_method) + 
    ggtitle(paste(current_method, "\nGreater in Controls"))  
  
  #Create merged graphs, assign and save graph 
  box_graph = ggarrange(fCD_graph, Control_graph, nrow = 2, common.legend = T)
  box_graph
  assign(paste(current_method, "graph", sep = "_"), box_graph)
}

#Arrange Graphs and save as jpeg 
plot_final = ggarrange(`Reads/g_graph`,`Relative Abundance_graph`,
          ncol =2 , common.legend = T); plot_final
jpeg ("PIME_RR_60.jpeg", res = 400, height = 4000, width = 4000)
plot_final
dev.off()
