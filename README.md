# Celiac_Disease_Pipeline

## Step 1: Subsetting Samples.  
  ***Input***: Phyloseq object, Metadata  
  ***Output***: Subset Metadata.csv, Subset Phyloseq object  
-  Overview:
   - Remove those with reads less than 1000 reads/g
   - Remove Non-celiac autoimmune and neurological  

## Step 2: Confounders 
  ***Input***: Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***Output***: Beta_Diversity.csv, Chisq_Pvalues.csv
-  Overview:
   - Find heterogeneous factors (Chisq.test)
   - Find Factors impacting Binomial Beta Diversity (Adonis)
   
## Step 3[.1]: [DESEQ](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  ***Input***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***Output***: Deseq_TotalAbundance_2.csv, Deseq_RelativeAbundance_2.csv.
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts

## Step 4: Merging Deseq output 
  ***Input***: Taxa Identifiers with 16S sequence, Deseq_TotalAbundance_2.csv<sub>(Step3)</sub>, Deseq_RelativeAbundance_2.csv<sub>(Step3)</sub>.  
  ***Output***: jpeg of Deseq LFC
-  Overview:
   - Remove Taxa found in fewer than the iteration threshold (50)
   - Graph LFC of all remaining iterations and taxa 
   
## Step 5: [PIME](https://github.com/microEcology/pime)
  ***Input***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***Output***: Deseq_TotalAbundance_2.csv, Deseq_RelativeAbundance_2.csv.
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts
   
## Step 3[.1]: [DESEQ](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  ***Input***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***Output***: Deseq_TotalAbundance_2.csv, Deseq_RelativeAbundance_2.csv.
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts
   
## Step 3[.1]: [DESEQ](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  ***Input***: Taxa Identifiers with 16S sequence, List of factors to match on <sub>(Step2)</sub>, Subset Phyloseq object <sub>(Step1)</sub>, Subset Sample data <sub>(Step1)</sub>.  
  ***Output***: Deseq_TotalAbundance_2.csv, Deseq_RelativeAbundance_2.csv.
-  Overview:
   - Run 100 iterations of matching (on the factors from step 2) and run through DESEQ
   - Total Abundance and Relative abundance in Different scripts
