# Celiac_Disease_Pipeline

## Step 1: Subsetting Samples.  
  ***INPUT***: Phyloseq object, Metadata  
  ***OUTPUT***: Subset Metadata.csv, Subset Phyloseq object  
-  Overview:
   - Remove those with reads less than 1000 reads/g
   - Remove Non-celiac autoimmune and neurological  

## Step 2: Confounders 
  ***INPUT***: Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***OUTPUT***: Beta_Diversity.csv, Chisq_Pvalues.csv
-  Overview:
   - Find heterogeneous factors (Chisq.test)
   - Find Factors impacting Binomial Beta Diversity (Adonis)
   
## Step 3[.1]: [DESEQ](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  ***INPUT***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***OUTPUT***: Deseq_TotalAbundance_2.csv, Deseq_RelativeAbundance_2.csv.
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts

## Step 4: Merging Deseq output 
  ***INPUT***: Taxa Identifiers with 16S sequence, Deseq_TotalAbundance_2.csv<sub>(Step3)</sub>, Deseq_RelativeAbundance_2.csv<sub>(Step3)</sub>.  
  ***OUTPUT***: jpeg of Deseq LFC
-  Overview:
   - Remove Taxa found in fewer than the iteration threshold (50)
   - Graph LFC of all remaining iterations and taxa 
   
## Step 5: [PIME](https://github.com/microEcology/pime)
  ***INPUT***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***OUTPUT***: PIME PCOA graphs, PIME taxa boxplots, 
-  Overview:
   - Begin by running 10 iterations to determine prevalence cutoff 
   - Run 100 iterations of matching and running through PIME, using the determined cutoff 
   - Create PCOA graphs of example iterations to show continual grouping: 25, 50, 75, 100
   - Create box plots of average reads for taxa in 50 or more iterations  
   
## Step 6: Prevalence 
  ***INPUT***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>
  ***OUTPUT***: Graph of prevalence (Absent in fCD, Absent in Controls, Large difference in prevalence)
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and find prevalence (/Percent of cohort with non-zero abundance)
   - Determine if absent in fCD, Controls, or difference > 25 
   - Graph if trend remains in 50 or more iterations 
   
## Step 7: Confounders
  ***INPUT***: Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>, List of important Deseq genera <sub>(Step4)</sub>, List of important PIME genera <sub>(Step6)</sub>
  ***OUTPUT***: CSVs (Genus, Variable, Pvalue, FDR), Graphs of DESEQ and PIME
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts
