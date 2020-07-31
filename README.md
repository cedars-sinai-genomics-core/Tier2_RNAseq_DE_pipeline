# get help manual
```bash
$ source activate RNAtier2
$ bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis.sh -h 
$ bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis_fixed_cutoff.sh -h
$ bash /home/genomics/genomics/apps/RNAseq_tier2/FFPE_RNAseq_analysis.sh -h
```
## Running pipeline using modified cutoff for comparisons
### 4 inputs are required for running RNAseq_analysis.sh
```bash
$ bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis.sh count.csv sample_info.csv comparison.csv project_ID
```
## Running pipeline using fixed cutoff for all comparisons
### 6 inputs are requried for running RNAseq_analysis_fixed_cutoff.sh
For example: 
####pvalue<0.05 was used as threshold
```bash
$ bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis_fixed_cutoff.sh count.csv sample_info.csv comparison.csv project_ID "pvalue<0.05" 1
```
####p_adj<0.01 and FC>=2 were used as threshold
```bash
$ bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis_fixed_cutoff.sh count.csv sample_info.csv comparison.csv project_ID "p_adj<0.01" 2
```
## for FFPE samples
```bash
$ bash /home/genomics/genomics/apps/RNAseq_tier2/FFPE_RNAseq_analysis.sh count.csv sample_info.csv comparison.csv project_ID
```


A few things to check before running the pipeline:  
1. The 1st column of sample_info file is colnames of COUNTS file.  
2. COUNTS file's column names should have same order as the 1st column of sample_info file.  
3. In sample_info file, Group_info (i.e. 3rd column) should be part of Sample_Name (i.e. 2nd column). Because Group_info will be used to grep what samples to include in the analysis.   
For example, if Group_info is G1 and Sample_Name should be G1_S1 or S1_G1. If Sample_Name is S1_Group1 or Group1_S1, no results will be generated.  

## Pathwat Analysis
```bash
/home/genomics/anaconda3/envs/RNAtier2/bin/Rscript ~/genomics/bin/Pathway_2020.R Corticotroph_KEGG.txt Corticotroph_GO_BP.txt Corticotroph_GO_CC.txt Corticotroph_GO_MF.txt Corticotroph
```
1. The 1st input is KEGG.txt. 
2. The next 3 inputs are GO terms txt (the order does not matter). 
3. The 5th input is SampleName (i.e. the outputs will be SampleName_DE_gene_KEGG.pdf, SampleName_DE_gene_GO_term.pdf, SampleName_DEGs_KEGG_enrichment.csv, SampleName_DEGs_GO_term_enrichment.csv). 
