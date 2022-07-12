
#####################
#Libraries and Functions 
#####################
library(vegan) #https://cran.r-project.org/web/packages/vegan/index.html
library(plyr) #http://myweb.facstaff.wwu.edu/minerb2/biometrics/plyr.html
library(phyloseq) #https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
library(ggplot2)
library(ggpubr)

'%notin%' = Negate('%in%')#https://stackoverflow.com/questions/38351820/negation-of-in-in-r

#####################
#Load Files 
#####################
map = read.csv("Celiac_Samples.csv")
PS = readRDS("Filters_PS.RDS")
taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")

#Create lists of Taxa of interests from steps 4 and 5
Deseq_list = c("Anaeroglobus", "Barnesiella","Candidatus_Soleaferrea", "Eubacterium", 
               "Monoglobus", "Senegalimassilia", "UCG.002")

PIME_list = c("Erysipelatoclostridium", "Haemophilus", "Lachnospiraceae_NK4A136_group", 
              "Butyricicoccus", "Collinsella", "Monoglobus", "Roseburia", "Ruminococcus", "Terrisporobacter")

################################
#Create list of columns to compare 
################################
environmental_list = c("DR13.DQ603","DR13.DQ604","DR14.DQ5", "DR14.DQ503","DR15.DQ601","DR15.DQ602","DR16.DQ5","DR16.DQ502","DR3.DQ2.5",
                       "DR4.DQ7","DR4.DQ8","DR5.DQ7","DR7.DQ2","DR7.DQ2.2","DR7.DQ2.5","DR7.DQ9","DR8.DQ4","DR9.DQ9","DQ2.5","DQ2.2","DQ8","DR1.DQ5",
                       "Antibiotics_during_Pregnancy","Infection_during_Pregnancy","Infection_as_newborn","Infection_w_antibiotics_1.12m",
                       "Cold.UpperRespTractInf","Otitis","Pneumonia","gastroenteritis","OtherInfection","OtherDisease",
                       "Antibiotics_Pregnancy","Corticosteroids_Pregnancy","HighBloodPressure_Meds_Pregnancy","Psychiatric_Meds_Pregnancy",
                       "PainKiller_Meds_Pregnancy","HormonePreparates_Pregnancy","Other_Meds_Pregnancy","Total_Meds_Pregnancy",
                       "Intro_CowsMilk_Binned","Intro_Gluten_Binned","Egg_first_year","Beef_first_year","Pork_first_year",
                       "Exclusive_Breastfeeding_Binned","Total_Breastfeeding_Binned","Intro_formula_Binned","Autoimmune_2_groups","Region","Sex","Siblings_at_birth","Delivery","Apartment","Father_Only_Elementary","Mother_Elementary",
                       "Father_Unemployed","Mother_Unemployed","Both_Parents_Abroad","Single_Mother","Stressful_Life_Event",
                       "Mother_NoSupport_Pregnancy", "Mother_NotSafe_Pregnancy","Worry_ChronicIllness_Child",
                       "Smoking_Pregnancy","Alcohol_Pregnancy_PerMonth","Risky_Alcohol_Pregnancy",
                       "Alcohol_Smoking_Medications_Pregnancy","Mother_Over_35","Father_Over_40" )

#####################
#Relative Abundance and Abundance 
ps.RA = transform_sample_counts(PS, function(x) x / sum(x) )

#################################################################
#qPCR; 16s rRNA reads/g
qpcr_data = map[,c("ID", "copies_16s_per_gram_stool")]
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
current_group = "PIME"

for (current_group in c("PIME", "DESEQ")) {
  final_data = data.frame()
  if (current_group == "PIME") {
    Genus_list = PIME_list
  } else {
    Genus_list = Deseq_list
  }
  
  for (current_method in c("Reads", "RelativeAbundace")) {
    #Read in appropriate Phyloseq object 
    if (current_method == "Reads") {
      current_ps = ps.RA.qpcr_1 
    } else {
      current_ps = ps.RA 
    }
    #####################
    #Create data.frame
    current_Taxa = data.frame(tax_table(current_ps))
    #Replace symbol with "." to save space 
    current_Taxa$Genus = gsub( "-","\\.", current_Taxa$Genus)
    current_Taxa = subset(current_Taxa, current_Taxa$Genus %in% Genus_list)
    tax_table(current_ps) = tax_table(as.matrix(current_Taxa))
    current_Taxa$ASV = rownames(current_Taxa)
    current_Taxa = current_Taxa[,c("ASV", "Genus")]
    
    #Create table of OTU counts and map data 
    count_table = data.frame(t(otu_table(current_ps)))
    count_table$ASV = row.names(count_table)
    count_table = merge(count_table, current_Taxa)
    row.names(count_table) = count_table[,"Genus"]
    count_table$ASV = NULL
    count_table[,"Genus"]= NULL
    count_table = data.frame(t(count_table))
    count_table$ID = row.names(count_table)
    count_table$ID = gsub("X", "", count_table$ID)
    count_table_2 = merge(count_table, map)
    
    if (current_method == "Reads") {
      Read_Count_table = count_table_2
    } else {
      RelAbun_Count_table = count_table_2
    }
    
    Output_data = data.frame()
    
    current_variable = "Sex"
    current_variable = environmental_list[1]
    Sign_data = data.frame()
    #Test significance of all metadata on Genera 
    for (current_variable in environmental_list) {
      pvalue_data = data.frame()
      current_data = count_table_2 
      current_data = subset(current_data, is.na(current_data[,current_variable]) == F)
      
      current_genus= Genus_list[1]
      #If there is more than one sub-category 
      if (length(unique(current_data[,current_variable]))>1) {
        for (current_genus in Genus_list) {
          #Run kruskal.test, extract and save p-value
          ANOVA = kruskal.test(as.formula(paste(current_genus, "~", current_variable)),
                               data = current_data)
          Pvalue = ANOVA$p.value
          P_row = data.frame("Genus" = current_genus, "Variable" = current_variable,  "Pvalue" =  Pvalue)
          pvalue_data = rbind(pvalue_data, P_row)
        }
      }
      #Adjust pvalues and subset 
      pvalue_data$FDR = p.adjust(pvalue_data$Pvalue, "fdr")
      pvalue_data = subset(pvalue_data, pvalue_data$Pvalue <= 0.05) 
      if (nrow(pvalue_data) > 0) {
        Sign_data = rbind(Sign_data, pvalue_data)
      }
    }
    Sign_data$method = current_method
    final_data = rbind(final_data, Sign_data)
  }
  
  #Create column of Genus and significant metadata variabke 
  final_data$Combined = paste(final_data$Genus, final_data$Variable, sep = "_")
  final_table = data.frame(table(final_data$Combined))
  #Keep those present in both relative abundance and total abundance 
  final_table = subset(final_table, final_table$Freq > 1)
  final_data = subset(final_data, final_data$Combined %in% final_table$Var1)
  
  Sign_data = final_data
  #####################################
  #Graphing  
  graphing_data = data.frame()
  graphing_data_1 = data.frame()
  current_variable = "Sex"
  current_genus = "Monoglobus"
  
  #Create long data for all significant variables 
  for (current_variable in unique(Sign_data$Variable)) {
    current_subset = subset(Sign_data,Sign_data$Variable == current_variable)
    for (current_genus in unique(current_subset$Genus)) {
      for (current_method in c("Reads", "RelAbun")) {
        if (current_method == "Reads") {
          current_count = Read_Count_table
        } else {
          current_count = RelAbun_Count_table
        }
        current_subset_1 = current_count[,c(current_variable, current_genus)]
        current_subset_1$Factor = current_variable
        current_subset_1$Genus = current_genus
        current_subset_1$Group = current_method
        colnames(current_subset_1) = c("Subfactor", "Value", "Factor", "Genus", "Group")
        current_subset_1= subset(current_subset_1, is.na(current_subset_1$Subfactor) == F)
        graphing_data = rbind(graphing_data, current_subset_1)
        current_subfactor = unique(current_subset_1$Subfactor)[1]
        for (current_subfactor in unique(current_subset_1$Subfactor)) {
          current_subset_3 = subset(current_subset_1, current_subset_1$Subfactor == current_subfactor )
          current_subset_3 = subset(current_subset_3, current_subset_3$Value > 0)
          current_mean = mean(current_subset_3$Value, na.rm = T)
          current_sd = sd(current_subset_3$Value, na.rm = T)
          current_row = data.frame(current_variable, current_genus, current_subfactor, 
                                   current_mean, current_sd, current_method)
          colnames(current_row) = c("Factor", 'Genus', "Subfactor", "Mean", "SD", "Group")
          graphing_data_1 = rbind(graphing_data_1, current_row)
        }
      }
    }
  }
  
  #Polish up column values 
  graphing_data$Factor = gsub("_", " ", graphing_data$Factor)
  graphing_data$Factor = ifelse(graphing_data$Factor == "Apartment", "Residence Type", graphing_data$Factor)

  graphing_data$Factor_Total = paste(graphing_data$Factor, graphing_data$Subfactor, sep = "\n")
  
  graphing_data$Subfactor = gsub(" ", "_", graphing_data$Subfactor)
  
  graphing_data$Genus = gsub("_group", " group", graphing_data$Genus)
  graphing_data$Genus = gsub("_", "\n", graphing_data$Genus)

  #Manually select colors for each subcatergory 
  break_list = c("Male","Female",
                 "No","Yes",
                 "1.2_times","3.5_times","Seldom", "Never",
                 "1.2_per_week","3.5_per_week","Daily",
                 "1_3","4_7","8_9",
                 "FALSE","TRUE",
                 "Other","Flat",
                 "0","1","2","3")
  value_list = c("cornflowerblue", "darksalmon", 
                 "tomato3", "springgreen4", 
                 "gray65", "gray55", "gray42", "gray22",
                 "plum2", "mediumorchid", "mediumorchid4", 
                 "deepskyblue", "dodgerblue", "dodgerblue4",
                 "red", "limegreen", 
                 "slateblue3", "turquoise3", 
                 "coral1", "brown2", "firebrick2", "firebrick4")

  #Create list of plots   
  plot_list = list()
  i = 1 #Create starting iteration value for plot list 
  for (current_genus in unique(graphing_data$Genus)) {
    current_subset = subset(graphing_data, graphing_data$Genus == current_genus)
    
    current_subset_RelAbun = subset(current_subset, current_subset$Group != "Reads")
    current_RelAbun_plot = ggplot(current_subset_RelAbun, aes(x = Factor_Total, y = Value, color = Subfactor)) +
      geom_boxplot() + 
      coord_cartesian(ylim = c(0, mean(current_subset_RelAbun$Value, na.rm = 2)*2)) +
      theme(legend.position="none",
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 5)) +
      scale_color_manual(breaks = break_list, values= value_list) +
      ggtitle(current_genus)
    current_plot = current_RelAbun_plot
    plot_list[[i]] = current_plot
    i = i + 1 
  }
  #Arrange all plots in list 
  final_plot = ggarrange(plotlist = plot_list, ncol = 3, nrow = 3)
  final_plot = annotate_figure(final_plot, paste(current_group, "\nRelative Abundance"))
  final_plot
  
  #Assign and save dataframe and plot
  assign(paste(current_group, "_plot", sep = ""), final_plot)
  print(paste(current_group, "_plot", sep = ""))
  assign(paste(current_group, "_data", sep = ""), Sign_data)
}

#Save dataframes of: Genus, Variable, Pvaluem FDR
write.csv(DESEQ_data, "DESEQ_confounder.csv", row.names = F)
write.csv(PIME_data, "PIME_confounder.csv", row.names = F)

#Save plots as Jpegs
jpeg("PIME_Confounders.jpeg", res = 400, height = 4000, width = 5000)
PIME_plot
dev.off()

jpeg("DESEQ_Confounders.jpeg", res = 400, height = 4000, width = 5000)
DESEQ_plot
dev.off()

