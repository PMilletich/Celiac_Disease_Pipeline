#PIME; https://github.com/microEcology/pime
library(phyloseq)
library(pime)
library(ggplot2)
library(ggpubr)
#library(vegan)

setwd("~/Desktop/ABIS/Celiac_Paper/Frontiers_Jun13/Reviewer_Edits_v2/")
Control_Num = 2

physeq_f = readRDS("Filters_PS.RDS") #[ 10389 taxa and 1756 samples ]
Sample_data = data.frame(sample_data(physeq_f))
taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
row.names(taxa_identifiers) = taxa_identifiers$ASV
taxa_identifiers = taxa_identifiers[,c("ASV", "Genus")]

#################################################################
#Relative Abundance 
#################################################################
ps.RA = transform_sample_counts(physeq_f, function(x) x / sum(x) )

#################################################################
#qPCR 
#################################################################
qpcr_data = Sample_data[,c("ID", "copies_16s_per_gram_stool")]
rownames(qpcr_data) = qpcr_data$ID
otu_RA = data.frame(otu_table(ps.RA))
otu_RA_1 = merge(qpcr_data, otu_RA, by = "row.names")
row.names(otu_RA_1) = otu_RA_1$Row.names
otu_RA_1$Row.names = NULL
otu_RA_1$ID = NULL
otu_RA_2 = otu_RA_1

for(i in 2:length(names(otu_RA_2))) {
  otu_RA_2[,i] <- round(otu_RA_2[,1] * otu_RA_2[, i])
}

otu_RA_2$copies_16s_per_gram_stool = NULL

ps.RA.qpcr = ps.RA
otu_RA_2[is.na(otu_RA_2)] = 0
otu_RA_2 = otu_RA_2[rowSums(otu_RA_2[])>0,]
otu_RA_2 = otu_RA_2[,colSums(otu_RA_2[])>0]

otu_table(ps.RA.qpcr) = otu_table(otu_RA_2, taxa_are_rows = FALSE)
ps.RA.qpcr

tax_table(ps.RA) = tax_table(ps.RA.qpcr)

seed_number = 1
current_method = "Reads"
subset_cases_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups != "Control")
subset_controls_OG = subset(Sample_data, Sample_data$Autoimmune_2_groups == "Control")

for (current_method in c("Reads", "RelativeAbundance")) {
  total_data = data.frame()
  best_df = data.frame()
  for (seed_number in 1:10) {
    current_run = paste(current_method, seed_number, sep = "_")
    iteration = seed_number
    seed_number = seed_number
    print(current_run)
    subset_cases = subset_cases_OG
    subset_controls = subset_controls_OG
    PIME_data = subset_cases
    PIME_data$Group = "fCD"
    
    Subset_1 = c("Region","Siblings_at_birth")
    current_id = "19289"#unique(subset_df_1$ID)[1]
    for (current_id in unique(subset_cases$ID)){
      current_row = subset(subset_cases, subset_cases$ID == current_id)
      current_row_subset = current_row[,Subset_1]
      current_row_subset = current_row_subset[ , colSums(is.na(current_row_subset)) == 0]
      
      if (current_id == "2061") {
        current_row_subset = data.frame("Region"= current_row_subset)
      }
      
      Matched_subset = merge(subset_controls, current_row_subset)
      if (nrow(Matched_subset) >= 1 & typeof(current_row_subset) == "list") {
        set.seed(seed_number)
        Matched_subset = Matched_subset[sample(nrow(Matched_subset), Control_Num), ]
        Matched_subset$Group = "Controls"
        PIME_data = rbind(PIME_data, Matched_subset)
        subset_controls = subset(subset_controls, ! subset_controls$ID %in% PIME_data$ID)
      } else {
        set.seed(seed_number)
        Matched_subset = subset_controls[sample(nrow(subset_controls), Control_Num), ]
        Matched_subset$Group = "Controls"
        PIME_data = rbind(PIME_data, Matched_subset)
        subset_controls = subset(subset_controls, ! subset_controls$ID %in% PIME_data$ID)
      }
    }
    table(PIME_data$Group)
    Sample_data = PIME_data
    rownames(Sample_data) = Sample_data$ID
    
    if (current_method == "Reads") {
      physeq_f = ps.RA.qpcr
    } else {
      physeq_f = ps.RA
    }
    sample_data(physeq_f) = sample_data(Sample_data)
    physeq_f
    ################################################
    #Celiac
    print(pime.oob.error(physeq_f, "Group"))
    per_variable_obj= pime.split.by.variable(physeq_f, "Group")
    
    print(sum(phyloseq::sample_sums(physeq_f)))
    
    prevalences=pime.prevalence(per_variable_obj)
    set.seed(42)
    best.prev=pime.best.prevalence(prevalences, "Group")
    best.df = data.frame(best.prev$`OOB error`)
    best.df$Iteration = seed_number
    if (seed_number == 1) {
      best_df= best.df
    } else {
      best_df = cbind(best_df, best.df)
    }
    imp60=best.prev$`Importance`$`Prevalence 60`
    prevalence.60 = prevalences$`60`
    input_ord = ordinate(prevalence.60, "PCoA" , "binomial")
    
    if (iteration%%25 == 0) {
      graph = plot_ordination(prevalence.60, input_ord , color = "Group")+
        stat_ellipse(aes(group = Group)) + 
        #facet_wrap(~Region) +
        ggtitle(paste("Iteration: ", iteration, sep = ""))
      assign(paste("graph_", iteration, sep = ""), graph)
    }
    
    imp60_subset = subset(imp60, imp60$MeanDecreaseAccuracy > 0)
    
    ###########
    #Data Fluffing
    count_table = data.frame(t(otu_table(prevalence.60)))
    count_table$ASV = row.names(count_table)
    count_table = merge(count_table, taxa_identifiers)
    row.names(count_table) = count_table$Genus
    Identifier_list = unique(rownames(count_table))
    tax_table_all_group = unique(count_table$Genus)
    count_table$ASV = NULL
    count_table$Genus= NULL
    count_table = data.frame(t(count_table))
    row.names(count_table) = gsub("X", "", row.names(count_table))
    count_table$ID = row.names(count_table)
    count_table_2 = merge(count_table, PIME_data)
    AI = subset(count_table_2, count_table_2$Autoimmune_2_groups != "Control")
    Control = subset(count_table_2, count_table_2$Autoimmune_2_groups == "Control")
    
    current_genus = unique(imp60_subset$Genus)[1]
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
    total_data = rbind(total_data, mean_data)
  }
  final_graph = ggarrange(graph_25, graph_50, 
                          graph_75, graph_100, 
                          ncol = 2, nrow = 2, common.legend = T)
  final_graph = annotate_figure(final_graph, current_method)
  assign(paste(current_method, "total", sep = ""), total_data)
  
  assign(paste("PCOA", current_method, sep = "_"), final_graph)
  
  
}
final_graph = ggarrange(`PCOA Reads _`, PCOA_RelativeAbundance)
jpeg(paste("./Images/", current_method, "_PIME_60.jpeg", sep = ""), res = 400, height = 2000, width = 4000)
final_graph
dev.off()

#write.csv(RelativeAbundancetotal, "./CSV_Output/RelAbund_PIME_100.csv", row.names = F)
#write.csv(Readstotal, "./CSV_Output/Reads_PIME_100.csv", row.names = F)

current_method = "Reads/g"
for (current_method in c("Reads/g", "Relative Abundance")) {
  if (current_method == "Reads/g") {
    total_data = Readstotal
  } else {
    total_data = RelativeAbundancetotal 
  }
  
  total_data_subset = subset(total_data, total_data$Controls >= 0 &
                               total_data$fCD > 0 & 
                               total_data$MeanDecreaseAccuracy > 0)
  
  count_data = data.frame(table(total_data_subset$Genus))
  count_data = subset(count_data, count_data$Freq >= 50)
  
  total_data_subset = subset(total_data_subset, total_data_subset$Genus %in% count_data$Var1)

  
  total_data_fCD = total_data_subset[,c("AI_Mean", "MeanDecreaseAccuracy", "Genus")]
  colnames (total_data_fCD) = c("Value", "Accuracy", "Genus")
  total_data_fCD$Group = "fCD"
  
  total_data_Control = total_data_subset[,c("Control_Mean", "MeanDecreaseAccuracy", "Genus")]
  colnames (total_data_Control) = c("Value", "Accuracy", "Genus")
  total_data_Control$Group = "Control"
  
  fCD_data = data.frame()
  Control_data = data.frame()
  current_genus = "Haemophilus"
  for (current_genus in unique(total_data_subset$Genus)) {
    AI_subset = subset(total_data_fCD, total_data_fCD$Genus == current_genus)
    AI_mean = mean(AI_subset$Value)
    
    Control_subset = subset(total_data_Control, total_data_Control$Genus == current_genus)
    Control_mean = mean(Control_subset$Value)
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
  
  box_graph = ggarrange(fCD_graph, Control_graph, nrow = 2, common.legend = T)
  box_graph
  assign(paste(current_method, "graph", sep = "_"), box_graph)
  
}
plot_final = ggarrange(`Reads/g_graph`,`Relative Abundance_graph`,
          ncol =2 , common.legend = T); plot_final

jpeg ("./Images/PIME_RR_60.jpeg", res = 400, height = 4000, width = 4000)
plot_final
dev.off()
