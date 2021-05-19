# This script is a stand-alone pathway analysis script that takes several command arguments
# you should first activate the RNAtier2 conda environment to make sure all of the necessary packages are loaded

### Example
# ./Pathway.R mouse DEG_list.csv Path/To/Results [LFC column number]



# Error Checking ----------------------------------------------------------
# Conda environment
if(system("which R", intern = T)!="/home/genomics/anaconda3/envs/RNAtier2/bin/R"){
  message("WARNING:\nYou should activate the RNAtier2 conda environment prior to running this script.")
}

## Command Args
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Speficy the input DEG list, species, and (optionally) a results path.", call.=FALSE)
}

# Set up ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse, quietly = T)
  library(clusterProfiler, quietly = T)
  library(magrittr, quietly = T)
  library(cowplot, quietly = T)
})

### Check species options
if(args[1]=="mouse"){
  suppressPackageStartupMessages(require(org.Mm.eg.db, quietly=T))
  message("Species set as ", args[1])
  assign("orgDB", org.Mm.eg.db)
  meta <- c(org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Mm.eg_dbInfo()[org.Mm.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "mmu")
} else if(args[1]=="human"){
  suppressPackageStartupMessages(require(org.Hs.eg.db, quietly=T))
  message("Species set as ", args[1])
  assign("orgDB", org.Hs.eg.db)
  meta <- c(org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Hs.eg_dbInfo()[org.Hs.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "hsa")
} else if(args[1]=="rat"){
  suppressPackageStartupMessages(require(org.Rn.eg.db, quietly = T))
  message("Species set as ", args[1])
  assign("orgDB", org.Rn.eg.db)
  meta <- c(org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Rn.eg_dbInfo()[org.Rn.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "rno")
} else if(args[1]=="dog"){
  suppressPackageStartupMessages(require(org.Cf.eg.db, quietly = T))
  message("Species set as ", args[1])
  assign("orgDB", org.Cf.eg.db)
  meta <- c(org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="ORGANISM",2], 
            org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="GOEGSOURCEDATE",2],            
            org.Cf.eg_dbInfo()[org.Cf.eg_dbInfo()[,1]=="KEGGSOURCEDATE",2],
            "cfa")
} else {
  stop("Species ", args[1], " not recognized. Specify human, mouse, rat, or dog.")
}

### Read in DEGs
tryCatch(DEG <- read_csv(args[2],
                         col_types = cols()),
         error=function(e){
           stop("Provide a valid CSV file of DEGs.")
         })
### Split combined ENSEMBL and Gene Symbol names if necessary
if(any(grepl(pattern = "_", x = DEG[1:5,]))){
  message("Parsing concatenated gene names into gene symbols")
  DEG %<>% 
    mutate(Row.names=str_split(Row.names, pattern = "_", simplify = T)[,2])
}
### Find log fold change column
if(!is.na(as.integer(args[4]))){
  message("Using user-specified column ", args[4]," as the log fold change column for GSEA.")
  LFC_index <- as.integer(args[4])
} else if(any(grepl(pattern = "log", colnames(DEG)))){
  LFC_index <- grep(pattern = "log", colnames(DEG))[1]
  LFC_name <- grep(pattern = "log", colnames(DEG), value = T)[1]
  message("Using column ", LFC_index, " ('",LFC_name,"') as the log fold change of expression for GSEA.\nIf this is incorrect, remove this column from the input csv file.")
} else {
  stop("Cannot find the log fold change of expression column for GSEA.\nCheck that it exists and/or provide a column number after the results path.")
}

message("Note: DEG list will be used as-is. No filtering of DEGs (e.g. by p-value) will be applied.")

### Set results path
if(str_sub(args[3],str_length(args[3]))!="/"){
  outdir <- paste0(args[3],"/")
} else {
  outdir <- args[3]
}
message("Outputting results to ",outdir)
dir.create(outdir,showWarnings = F)

# Run GSEA analysis ------------------------------------------------------------
s1 <- DEG %>% 
  dplyr::arrange(-.[,LFC_index]) %>% 
  dplyr::select(all_of(c(1,LFC_index))) %>% 
  deframe() %>% 
    gseGO(ont = "ALL",
          OrgDb = orgDB,
          keyType = "SYMBOL",
          pvalueCutoff = Inf,
          eps = 1e-200) %>% 
  data.frame()

# Write CSV ---------------------------------------------------------------
if(dim(s1)[1]>0){
  s2 <- s1 %>% 
    dplyr::arrange(pvalue) %>% 
    dplyr::select(ONTOLOGY, ID, Description, NES, setSize, pvalue, p.adjust, core_enrichment) %T>% 
    write_csv(paste0(outdir,"GSEA.csv"))
} else {
  message("No enriched GO terms found.\nNo GO plots or lists will be generated.")
}

# Plot --------------------------------------------------------------------
if(exists("s2")){
  # raw pvalue
  s3 <- s2 %>% 
    dplyr::filter(pvalue<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(.$NES)]))
  if(dim(s3)[1]>0){
    s3 %>% ggplot(aes_string(x="NES",
                             y="Concatenated",
                             size="setSize",
                             colour="pvalue"))+
      facet_grid(ONTOLOGY~.,
                 scales="free_y",space="free_y")+
      geom_point()+
      scale_colour_gradient(low="red",
                            high="blue",
                            name="Raw p-value")+
      scale_size(name="Num. Genes")+
      labs(x="Normalized Enrichment Score",
           y="Term",
           title="Enriched GO Terms",
           caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[2]))) +
      theme_minimal_grid() +
      theme(strip.background = element_rect(
        color="black",
        fill="#CCCCCC",
        size=1,
        linetype="solid"))
    ggsave2(filename = paste0(outdir,"GSEA_p0.05.pdf"),
            height=dim(s3)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else {
    message("No significant GO terms found at p<0.05.\nNo plots will be generated.")
    dev.off()
  }
  # adjusted pvalue
  s4 <- s2 %>% 
    dplyr::filter(p.adjust<0.05) %>% 
    mutate(Concatenated = paste0(ID,"~",Description)) %>% 
    mutate(Concatenated = str_trunc(Concatenated, 60)) %>% 
    mutate(Concatenated = factor(Concatenated,
                                 levels = .$Concatenated[order(.$NES)]))
  if(dim(s4)[1]>0){
    s4 %>% 
      ggplot(.,aes_string(x="NES",
                           y="Concatenated",
                           size="setSize",
                           colour="p.adjust"))+
          facet_grid(ONTOLOGY~.,
                     scales="free_y",space="free_y")+
          geom_point()+
          scale_colour_gradient(low="red",
                                high="blue",
                                name="Adj. p-value")+
          scale_size(name="Num. Genes")+
          labs(x="Normalized Enrichment Score",
               y="Term",
               title="Enriched GO Terms",
               caption=bquote(italic(.(meta[1]))*', Sourced:'~.(meta[2]))) +
      theme_minimal_grid() +
          theme(strip.background = element_rect(
            color="black",
            fill="#CCCCCC",
            size=1,
            linetype="solid"))
    ggsave2(filename = paste0(outdir,"GSEA_padj0.05.pdf"),
            height=dim(s4)[1]*0.25+2.5,
            width=10,
            limitsize = F)
  } else
    message("No significant GO terms found at adjusted p<0.05.\nNo plot will be generated.")
}

# No idea why this is create, but this will remove it
if(file.exists("./Rplots.pdf")){
  file.remove("./Rplots.pdf")
}

