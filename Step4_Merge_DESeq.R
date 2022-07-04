library(ggplot2)
library(ggpubr)

taxa_identifiers = read.csv("Taxa_ASV_Identifiers.csv")
taxa_identifiers_1 = data.frame()
current_glom = "Phylum"
for (current_glom in c("Genus")) {
  current_subset = unique(taxa_identifiers[,current_glom])
  current_df = data.frame("Taxa" = current_glom, "Taxa_Name" = current_subset)
  taxa_identifiers_1= rbind(taxa_identifiers_1, current_df)
}

taxa_identifiers_1$Taxa_Name = gsub("-", ".", taxa_identifiers_1$Taxa_Name)
taxa_identifiers_1$Taxa_Name = gsub("/", ".", taxa_identifiers_1$Taxa_Name)

iteration_Threshold = 50
current_method = "Reads/g"

for (current_method in c("Reads/g", "Relative Abundance")) {
  if (current_method == "Reads/g") {
    current_data = read.csv("./CSV_Output/Deseq_Matched_2.csv")
  } else {
    current_data = read.csv("./CSV_Output/Deseq_MatchedRelAbun_2.csv")
  }
  
  df_Count = data.frame(table(current_data$Taxa_Name))
  df_Count = subset(df_Count, df_Count$Freq >= iteration_Threshold)
  
  Genus_list = unique(df_Count$Var1)
  
  current_data = subset(current_data, current_data$Taxa_Name %in% Genus_list)
  df_Count = data.frame(table(current_data$Taxa_Name))
  colnames(df_Count) = c("Taxa_Name", "Iterations")
  current_data_Condensed = data.frame()
  
  current_data = merge(df_Count, current_data)
  
  current_data$Taxa = NULL
  
  current_Taxa = unique(current_data$Taxa_Name)[1]
  for (current_Taxa in unique(current_data$Taxa_Name)) {
    current_subset = subset(current_data, current_data$Taxa_Name == current_Taxa)
    
    current_LFC = mean(current_subset$log2FoldChange)
    current_LFC_SD = sd(current_subset$log2FoldChange)
    
    current_AI_prev = mean(current_subset$AI_prev)
    current_AI_prev_SD = sd(current_subset$AI_prev)
    current_Control_prev = mean(current_subset$Control_prev)
    current_Control_prev_SD = sd(current_subset$Control_prev)
    
    current_row = data.frame("Taxa_Name" = unique(current_subset$Taxa_Name),
                             "Iteration" = unique(current_subset$Iterations),
                             "LFC" = current_LFC, "LFC_SD" = current_LFC_SD,
                             "AI_Prev" = current_AI_prev, "AI_Prev_SD" = current_AI_prev_SD, 
                             "Control_Prev" = current_Control_prev, "Control_Prev_SD" = current_Control_prev_SD,
                             "Group" = current_method)
    
    current_data_Condensed = rbind(current_data_Condensed, current_row)
  }
  
  current_data = current_data_Condensed
  current_data$Taxa_Name = gsub("_", " ", current_data$Taxa_Name)
  long_data = current_data
  
  long_data$Trend = ifelse(long_data$LFC < 0, "Controls", 
                           ifelse(long_data$LFC > 0, "fCD", "NA"))
  
  #long_data = long_data[complete.cases(long_data),]
  print(current_method)
  print(unique(long_data$Taxa_Name))
  LFC = ggplot(long_data, aes(x = LFC, y = Taxa_Name)) + 
    geom_point( position = position_dodge(width=0.9)) + 
    geom_errorbar(position = position_dodge(width=0.9),
                  aes(y=Taxa_Name, xmin=LFC-LFC_SD, xmax=LFC+LFC_SD), 
                  width=0.4, alpha=0.9, size=1.3) +
    geom_vline(xintercept = 0) +
    theme(axis.title.y = element_blank()) +
    ggtitle(paste("Log2FoldChange", current_method, sep = "\n"))
  LFC
  
  if (current_method == "Reads/g") {
    current_method = "Reads"
  } else {
    current_method= "RelAbun"
  }
  jpeg(paste("./Images/Deseq_50_", current_method, ".jpeg", sep = ""), res = 400, height= 3000, width = 3000)
  print(LFC)
  dev.off()
  assign(paste(current_method, "_plot", sep = ""),LFC )
  print(paste(current_method, "_plot", sep = ""))
}

jpeg("./Images/Deseq_Combined.jpeg", res = 400, height = 2000, width = 4000)
ggarrange(Reads_plot, RelAbun_plot, ncol = 2)
dev.off()
