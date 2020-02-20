#########################################################################
# File Name: RNAseq_analysis_fixed_cutoff.sh
# Author: Di Wu
# mail: di.wu@cshs.org
# Created Time: Feb., 2020
#########################################################################
#!/bin/bash
#!/bin/bash
display_usage() {
	echo -e "NAME:\n  RNAseq_analysis."
	echo -e "\nDESCRIPTION:\n   This pipeline will use COUNT file to do doenstream differentially expressed genes (DEGs) analysis."
	echo -e "\nUsage:\n   bash /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_analysis.sh Count_file.csv Sample_info.csv comparisons.csv Ptoject_ID "pvalue<0.05" 1"

    echo "Input options:"
    echo "   -h|--help    show this help"

    echo "Input files:"
    echo "   The first input is Count_file.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_COUNTS.csv  as an example format"
    echo "   The second input is Sample_info.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_sample_info.csv as an example format"
    echo "   The third input is Conparisons.csv, check https://github.com/dxw5099/RNAseq_downstream/blob/master/demo_comparisons.csv as an example format"
    echo "   The fourth input is project_ID, for example: AP-5782–11–08–2018"
    echo "   The fifth input is p-value or adjusted p-value, for example: p_adj<0.01 or pvalue<0.05"
    echo "   The sixth input is FC value, for example: 1 or 2 or 3..."
	  echo "                                                  1 means no filtering"
    echo "                                                  2 means FoldChange >=2"
    echo "                                                  3 means FoldChange >=3"
    }
# check whether user had supplied -h or --help . If yes display usage
if [[ ($1 == "--help") ||  ($1 == "-h") ]]
then
	display_usage
	exit 0
fi

# get PCA plots for all samples and DEGs table, interactive report for each comparison
#/hpc/apps/R/3.4.1/bin/Rscript /common/genomics-core/data/Temp/Di_RNA_seq_test/downstream_test/RNAseq_tier2.R $1 $2 $3 $4
/usr/bin/Rscript /home/genomics/genomics/apps/RNAseq_tier2/RNAseq_tier2_fixed_cutoff.R $1 $2 $3 $4 $5 $6
