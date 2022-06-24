library(phyloseq)

setwd("~/Desktop/Reviewer_Edits_V2/")

'%notin%' = Negate('%in%')

#Methods
OG_PS = readRDS("phyloseq.13.Final.RDS") 
Start = read.csv("Regions_MetaData_2022_Neuro.csv")
rownames(Start) = Start$ID  
sample_data(OG_PS) = sample_data(Start)

otu_21 = data.frame(otu_table(OG_PS))
otu_21$Counts = rowSums(otu_21)
summary(otu_21$Counts)
otu_212 = subset(otu_21, otu_21$Counts > 1000)
otu_table(OG_PS) = otu_table(otu_212, taxa_are_rows = F)

With_Gut_Celiac_only_1000 = data.frame(sample_data(OG_PS))
With_Gut_Celiac_only_1000 = subset(With_Gut_Celiac_only_1000, is.na(With_Gut_Celiac_only_1000$copies_16s_per_gram_stool) == F)

OG_PS
sample_data(OG_PS) = sample_data(With_Gut_Celiac_only_1000)
OG_PS #1346 Samples

table(With_Gut_Celiac_only_1000$Autoimmune_Disorder)

###########################################################################
#Remove Autoimmune 
###########################################################################
Autoimmune = subset(With_Gut_Celiac_only_1000, With_Gut_Celiac_only_1000$Autoimmune_Disorder != "Control")
Autoimmune_NonCeliac = subset(Autoimmune, grepl("Celiac", Autoimmune$Autoimmune_Disorder) == F)
nrow(Autoimmune_NonCeliac)

Autoimmune_Celiac = subset(Autoimmune, grepl("Celiac", Autoimmune$Autoimmune_Disorder) == T)

With_Gut_Celiac_only = subset(With_Gut_Celiac_only_1000, grepl("Celiac", With_Gut_Celiac_only_1000$Autoimmune_Disorder) == T |
                                With_Gut_Celiac_only_1000$Autoimmune_Disorder == "Control")

rownames(With_Gut_Celiac_only) = With_Gut_Celiac_only$ID
sample_data(OG_PS) = sample_data(With_Gut_Celiac_only)

table(With_Gut_Celiac_only$Autoimmune_Disorder)


###########################################################################
#Remove Neurological  
###########################################################################
Neurological = subset(With_Gut_Celiac_only, With_Gut_Celiac_only$Neurological == "Control" | 
                        With_Gut_Celiac_only$Autoimmune_Disorder == "Celiac")

rownames(Neurological) = Neurological$ID
sample_data(OG_PS) = sample_data(Neurological)

table(Neurological$Autoimmune_Disorder)

ps.filtered = filter_taxa(OG_PS, function(x) sum(x > 5) >= 5, TRUE)
ps.filtered = tax_glom(ps.filtered, "Genus")

otu_RA_2 = data.frame(otu_table(ps.filtered))
current_ASV = colnames(otu_RA_2)[1]
# for (current_ASV in colnames(otu_RA_2)) {
#   current_list = otu_RA_2[,current_ASV]
#   Prevalence = (length(subset(current_list, current_list > 0))/length(current_list))*100
#   if (Prevalence < 1) {
#     otu_RA_2[,current_ASV] = NULL
#   }
# }
# dim(otu_RA_2)

otu_table(ps.filtered) = otu_table(otu_RA_2, taxa_are_rows = F)


otu_21 = data.frame(otu_table(ps.filtered))
otu_21$Counts = rowSums(otu_21)
summary(otu_21$Counts)
otu_212 = subset(otu_21, otu_21$Counts > 1000)


write.csv(Neurological, "Celiac_Samples.csv", row.names = F)
saveRDS(ps.filtered, "Filters_PS.RDS")
