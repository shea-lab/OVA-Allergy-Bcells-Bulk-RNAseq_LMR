

# Code Written by Laila M Rad
# Lonnie Shea Lab
# Updated July 2025


#load libraries ----
pacman::p_load(tidyverse, plyr, magrittr, stats, dplyr, limma, RColorBrewer, gplots, glmnet, 
               biomaRt, colorspace, ggplot2, fmsb, car, mixOmics, DESeq2, apeglm, rgl, qpcR,
               boot, caret, ggvenn, grid, devtools, reshape2, gridExtra, factoextra, edgeR, 
               cowplot, pheatmap, coefplot, randomForest, ROCR, genefilter, Hmisc, rdist, 
               factoextra, ggforce, NormqPCR, ggpubr, matrixStats, GSEAmining, ggrepel)

library(ggbiplot)
library(ggprism)
library(tibble)
library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(msigdbr)
library(pathview)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(enrichR)
library(EnhancedVolcano)
library("DEGreport")
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GOSemSim)
library(ggtreeExtra)
library(ggtree)
#BiocManager::install("ggtreeExtra")


# Read in data ----
# Read in the file. Don't set row names yet
# Note if using R < 4.0.0, set stringsAsFactors = FALSE in read.delim
#dataset location
path = '/Users/lailarad/Documents/UMich Research/Allergy/cohort 14 APC:Bcell co culture/B cell RNA seq/' 

#set working directory to dataset location
setwd(path)
getwd()
#find dataset in location
file = list.files(path=path,
                  pattern = ".txt" )
file

#read in count data ----
data <- read.delim("gene_expected_count.annot.txt", 
                   row.names = NULL)


#data <- read.delim("deliverables/counts/gene_expected_count.annot.txt", row.names = NULL)
# Deal with genes that don't have annotated gene symbols (external_gene_name)
# Use ENSEMBL ID if gene symbol not available
data$external_gene_name <- ifelse(
  data$external_gene_name == ".",
  data$gene_id,
  data$external_gene_name
)

# Deal with duplicated gene symbols
# Combine gene symbol with ENSEMBL ID if non-unique
data$external_gene_name <- ifelse(
  duplicated(data$external_gene_name),
  paste(data$external_gene_name, data$gene_id, sep="_"),
  data$external_gene_name
)

# Then we can use the gene symbol column as the row names,
# and subset the count data for further analysis
rownames(data) <- data$external_gene_name
count.data <- data[,5:ncol(data)] # All columns after 4 are count data
gene_info = data[,1:4] # Columns with gene names and description 

#renumber sample 19 to be 11 (resubmitted sample 11 with sample 19)
colnames(count.data) <- sub("X8188.LR.19", "X8188.LR.11", colnames(count.data))

#just leave number in sample name 
colnames(count.data) <- sub("X8188.LR.", "", colnames(count.data))

#sample information
#read information about samples provided by AGS Core
info_file1= read.csv('DemuxStats_8188-LR.csv', 
                     check.names=FALSE,header=T, sep=",")

#remove resubmit from sample name 
info_file1$Description <- sub("_resubmit", "", info_file1$Description)

#read in sample information about experimental conditions
info_file2 = read.csv('cohort14_bcells_rnaseq_sample_labels.csv', 
                      check.names=FALSE,header=T, sep=",")

#make dataframe with sample info
sample_info = merge(info_file1, info_file2, by = "Description")

#remove sample name to leave only sample ID description
sample_info$SampleNum = sub("FA_cohort14_sID_", "", sample_info$Description) 

sample_info[which(sample_info$Tx == "PLG(OVA-hi) NP"),"Tx"] = "OVA NPs"

sample_info$Tissue_Tx = paste0(sample_info$Tissue, "_", sample_info$Tx)

rownames(sample_info) <- sample_info$SampleNum
#write.csv(sample_info, "sample_info_fromR.csv")


#check the order information
count.data <- count.data[, rownames(sample_info)] #reorder to match
all(rownames(sample_info) == colnames(count.data)) #Needs to be true to be ordered correctly

#filtering ----
filt_low <- 10 # These three filters are per sample per gene

#Spleen count data for PBS, OVA NP, Soluble OVA treated mice
Spleen = count.data[,c(sample_info$Tissue == "Spleen") ]

#meta data for spleen data
Spleen_info = sample_info[which(sample_info$Tissue == "Spleen"  ),]

#PP count data for PBS and OVA NP treated mice
PP_NP = count.data[,c(sample_info$Tissue == "PP" & sample_info$Tx %in% c("PBS", "OVA NPs"))]
PP_sOVA = count.data[,c(sample_info$Tissue == "PP" & sample_info$Tx %in% c("PBS", "Soluble OVA"))]
PP_NPvsOVA = count.data[,c(sample_info$Tissue == "PP" & sample_info$Tx %in% c("OVA NPs", "Soluble OVA"))]

#meta data for PP data
PP_NP_info = sample_info[which(sample_info$Tissue == "PP" & sample_info$Tx %in% c("PBS", "OVA NPs")),]
PP_sOVA_info = sample_info[which(sample_info$Tissue == "PP" & sample_info$Tx %in% c("PBS", "Soluble OVA")),]
PP_NPvsOVA_info = sample_info[which(sample_info$Tissue == "PP" & sample_info$Tx %in% c("OVA NPs", "Soluble OVA")),]

sample_num = length(Spleen) #number of spleen samples
Spleen_filt <- Spleen[rowSums(Spleen)>0,] # Remove 0 expression genes

#keep rows that have less than sample_num*(3/4) [6*3/4 = 4.5] number of zero counts in the row
Spleen_filt <- Spleen_filt[rowSums(Spleen_filt == 0) <= sample_num*(3/4),] # Remove mostly 0 genes

# Remove low expression genes
keep1 <- rowSums(Spleen_filt) > (filt_low*sample_num) # Low count filter, required
Spleen_filt <- Spleen_filt[keep1,] # Apply filter to counts data
dim(Spleen_filt) # Dimensions after applying low filter

is.na(Spleen_filt) %>% table() # Check for NA values

sample_num = length(PP_NP) #number of PP samples
PP_NP_filt <- PP_NP[rowSums(PP_NP)>0,] # Remove 0 expression genes

#keep rows that have less than sample_num*(3/4) [6*3/4 = 4.5] number of zero counts in the row
PP_NP_filt <- PP_NP_filt[rowSums(PP_NP_filt == 0) <= sample_num*(3/4),] # Remove mostly 0 genes

# Remove low expression genes
keep1 <- rowSums(PP_NP_filt) > (filt_low*sample_num) # Low count filter, required
PP_NP_filt <- PP_NP_filt[keep1,] # Apply filter to counts
dim(PP_NP_filt) # Dimensions after applying low filter
is.na(PP_NP_filt) %>% table() # Check for NA values

sample_num = length(PP_sOVA) #number of PP samples
PP_sOVA_filt <- PP_sOVA[rowSums(PP_sOVA)>0,] # Remove 0 expression genes

#keep rows that have less than sample_num*(3/4) [6*3/4 = 4.5] number of zero counts in the row
PP_sOVA_filt <- PP_sOVA_filt[rowSums(PP_sOVA_filt == 0) <= sample_num*(3/4),] # Remove mostly 0 genes

# Remove low expression genes
keep1 <- rowSums(PP_sOVA_filt) > (filt_low*sample_num) # Low count filter, required
PP_sOVA_filt <- PP_sOVA_filt[keep1,] # Apply filter to counts
dim(PP_sOVA_filt) # Dimensions after applying low filter
is.na(PP_sOVA_filt) %>% table() # Check for NA values

sample_num = length(PP_NPvsOVA) #number of PP samples
PP_NPvsOVA_filt <- PP_NPvsOVA[rowSums(PP_NPvsOVA)>0,] # Remove 0 expression genes

#keep rows that have less than sample_num*(3/4) [6*3/4 = 4.5] number of zero counts in the row
PP_NPvsOVA_filt <- PP_NPvsOVA_filt[rowSums(PP_NPvsOVA_filt == 0) <= sample_num*(3/4),] # Remove mostly 0 genes

# Remove low expression genes
keep1 <- rowSums(PP_NPvsOVA_filt) > (filt_low*sample_num) # Low count filter, required
PP_NPvsOVA_filt <- PP_NPvsOVA_filt[keep1,] # Apply filter to counts
dim(PP_NPvsOVA_filt) # Dimensions after applying low filter
is.na(PP_NPvsOVA_filt) %>% table() # Check for NA values


#DESeq2 ----
#make Treatment conditions as factors
Spleen_info$Tx = factor(Spleen_info$Tx, levels = c("PBS", "OVA NPs", "Soluble OVA"))
PP_NP_info$Tx = factor(PP_NP_info$Tx, levels = c("PBS", "OVA NPs"))

#Remove genes Gm- or -Rik or -ps
Spleen_filt1 = Spleen_filt[!grepl("^Gm|Rik$|Rik[0-9]+$|-ps[0-9]+$|-ps$", rownames(Spleen_filt)),] 
PP_NP_filt1 = PP_NP_filt[!grepl("^Gm|Rik$|Rik[0-9]+$|-ps[0-9]+$|-ps$", rownames(PP_NP_filt)),] 

dds_full_Spleen = DESeqDataSetFromMatrix(
  countData = Spleen_filt1,
  colData = Spleen_info,
  design = ~ Tx
)
dds_full_Spleen$Tx
dds_full_Spleen$Tissue
dds_full_Spleen = DESeq(dds_full_Spleen)

resultsNames(dds_full_Spleen)

dds_NP_PP = DESeqDataSetFromMatrix(
  countData = PP_NP_filt1,
  colData = PP_NP_info,
  design = ~ Tx
)
dds_NP_PP$Tx
dds_NP_PP$Tissue
dds_NP_PP = DESeq(dds_NP_PP)

resultsNames(dds_NP_PP)




LFc = 0.25 #logfold change cutoff for DEGs
pc = 0.1 #pvalue cutoff  for DEGs




#DEGs in Spleen----
#Spleen_NPvPBS = results(dds_full_Spleen, contrast = c("Tx", "OVA NPs", "PBS"), alpha = pc )
#Spleen_NPvPBS
#Log fold change shrinkage for visualization and ranking----
#shrinkage methods provided by DESeq2 are good for ranking genes by “effect size”
#Spleen_NPvPBS1 = lfcShrink(dds_full_Spleen, coef = "Tx_OVA.NPs_vs_PBS", type ="apeglm")

#Spleen_NPvPBS2 = lfcShrink(dds_full_Spleen, res = Spleen_NPvPBS, type="ashr")

#an adaptive Normal distribution as prior shrinkage estimator 
Spleen_NPvPBS3 = lfcShrink(dds_full_Spleen, contrast = c("Tx", "OVA NPs", "PBS"), type="normal")
Spleen_NPvPBS3
summary(Spleen_NPvPBS3)

genesofinterest = c("Fcer2a", "Elane", "Mcpt1", "Mcpt2", "Itgb1", "Cd19")
Spleen_NPvPBS3[genesofinterest,]
#Spleen_NPvPBS2[genesofinterest,]
#Spleen_NPvPBS1[genesofinterest,]


#Upreg genes in OVA NP in Spleen
SpleenUpNP = Spleen_NPvPBS3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < pc, log2FoldChange > LFc) %>%
  arrange(desc(log2FoldChange))

# Upreg genes in PBS in Spleen (Downreg in OVA NPs in Spleen )
SpleenDownNP = Spleen_NPvPBS3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < pc, log2FoldChange < -LFc) %>%
  arrange(desc(abs(log2FoldChange)))


#DEGs in PP----
#Log fold change shrinkage for visualization and ranking----
#shrinkage methods provided by DESeq2 are good for ranking genes by “effect size”

#PP_NPvPBS = results(dds_NP_PP, contrast = c("Tx", "OVA NPs", "PBS"), alpha = 0.1, pAdjustMethod = "BH")
#summary(PP_NPvPBS)
#deseq_res <- ihw(PP_NPvPBS, alpha=0.1, adjustment_type = "BH") #ihw method
#summary(deseq_res)
#PP_NPvPBS3ihw = lfcShrink(dds_NP_PP, contrast = c("Tx", "OVA NPs", "PBS"), type="normal", res = deseq_res)
#summary(PP_NPvPBS3ihw)

PP_NPvPBS3 = lfcShrink(dds_NP_PP, contrast = c("Tx", "OVA NPs", "PBS"), type="normal")
summary(PP_NPvPBS3)


#Upreg genes in OVA NP in PP
PPUpNP = PP_NPvPBS3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < pc, log2FoldChange > LFc) %>%
  arrange(desc(abs(log2FoldChange)))

#Upreg genes in PBS in PP (Down in OVA NP)
PPDownNP = PP_NPvPBS3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < pc, log2FoldChange < -LFc) %>%
  arrange(desc(abs(log2FoldChange)))


#Lit review Breg markers ----
Bcellmarkers = read.csv("Bcellmarkers2.csv")

#Genes upregulated in Bregs
BregUp = Bcellmarkers[which(Bcellmarkers$Regulation == "Up" & Bcellmarkers$Subset == "Regulatory"),]
#Genes downregulated in Bregs
BregDown = Bcellmarkers[which(Bcellmarkers$Regulation == "Down" & Bcellmarkers$Subset == "Regulatory"),]


intersect(BregUp$Genes, SpleenUpNP$gene)
intersect(BregUp$Genes, SpleenDownNP$gene)
intersect(BregDown$Genes, SpleenDownNP$gene)
intersect(BregDown$Genes, SpleenUpNP$gene)


intersect(BregUp$Genes, SpleenUpNP$gene)




intersect(BregDown$Genes, SpleenDownNP$gene)



intersect(BregUp$Genes, PPUpNP$gene)
intersect(BregUp$Genes, PPDownNP$gene)
intersect(BregDown$Genes, PPDownNP$gene)
intersect(BregDown$Genes, PPUpNP$gene)

#Breg Markers in DEGs ----

# DESeq2 'regularized log' transformation
# transforms the count data to the log2 scale in a way which minimizes 
# differences between samples for rows with small counts, 
# and which normalizes with respect to library size. 
rld_Sp <- rlog(dds_full_Spleen, blind=FALSE)
rld_PP <- rlog(dds_NP_PP, blind=FALSE)

# the matrix of transformed values is stored in assay(rld)
rld_mat_Sp = assay(rld_Sp)
rld_mat_PP = assay(rld_PP)

plotPCA(rld_Sp, intgroup = "Tx")
plotPCA(rld_PP, intgroup = "Tx")

# Obtain rlog values for significant DEGs
cluster_rlog <- rld_mat_Sp[c(SpleenUpNP$gene, SpleenDownNP$gene), ]

# Breg markers that are DEGs
BregGenesMap = c(
  intersect(BregUp$Genes, SpleenUpNP$gene),
  intersect(BregUp$Genes, SpleenDownNP$gene),
  intersect(BregDown$Genes, SpleenDownNP$gene)
)

BregMarkersDEGsSp = intersect(
    #DEGs in Spleen
    Spleen_NPvPBS3 %>%
      data.frame() %>%
      rownames_to_column(var="gene") %>% 
      as_tibble() %>% 
      filter(padj < pc, abs(log2FoldChange) > LFc) %>%
      dplyr::select(gene) %>%
      deframe(),
      #Breg Markers
      Bcellmarkers[which( Bcellmarkers$Subset == "Regulatory"),"Genes"]
)


# Remove duplicate genes
BregGenesMap = unique(BregGenesMap)

# Save information about Breg Markers that are DEGs in a seperate csv file
#write.csv( Bcellmarkers[which(Bcellmarkers$Genes %in% BregGenesMap & Bcellmarkers$Subset == "Regulatory"), ], 
#           "BregMarkersDEGsfromR4.csv")

# Dataframe with only the Breg markers that are DEGs
BregmarkersDEG = read.csv("BregMarkersDEGsfromR4.csv")
head(BregmarkersDEG)

colnames(BregmarkersDEG)

setdiff(BregmarkersDEG$Genes, BregGenesMap)
setdiff(BregGenesMap, BregmarkersDEG$Genes)

setdiff(BregMarkersDEGsSp, BregmarkersDEG$Genes)
setdiff(BregmarkersDEG$Genes, BregMarkersDEGsSp)

#**plot heatmap of Breg Marker DEGs ----

# Degs 0.25 log2 fc and padj less than 0.1;
# Heatmap of log-transformed values across samples
# row scaled by removing the mean (centering) and dividing by the standard deviation (scaling) per tissue

# Information about samples for plotting saved in df dataframe
df1 <- as.data.frame(colData(dds_full_Spleen)[,c("Tissue","Tx")])

rownames(df1) = colnames(cluster_rlog)

df1 = df1 %>%
  filter(Tissue == "Spleen", Tx %in% c("OVA NPs", "PBS"))

df1

cluster_rlog_Sp = as.data.frame(cluster_rlog[, rownames(df1)])

cluster_rlog_Sp$gene = rownames(cluster_rlog_Sp)

cluster_rlog_Sp <- inner_join(cluster_rlog_Sp, BregmarkersDEG[BregmarkersDEG$Genes %in% unique(c(intersect(BregUp$Genes, SpleenUpNP$gene),
                                                                                                 intersect(BregUp$Genes, SpleenDownNP$gene),
                                                                                                 intersect(BregDown$Genes, SpleenDownNP$gene))
                                                                                               ),], by=c("gene"="Genes")) 
rownames(cluster_rlog_Sp) = cluster_rlog_Sp$gene
head(cluster_rlog_Sp)

cluster_rlog_Sp = cluster_rlog_Sp %>%
  arrange(Regulation)
colnames(cluster_rlog_Sp)[10] = "Expression in Bregs"

# colors for plotting
anno_colors = list(
  #Tissue = c(Spleen ="orange1", PP = "purple3"),
  Tx = c("PBS" = "grey", "OVA NPs" = "#89CFF0"),
  #Regulation = c(Up ="firebrick3", Down = "blue2", `Up or Down` = "purple2")
  `Expression in Bregs` = c(Up ="firebrick3", Down = "blue2", `Up or Down` = "purple2")
)





# Heatmap of Breg Markers DEGs in Spleen ----

#png(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.png",units = "in",width = 3.5, height = 6.2, res = 400)
#cairo_pdf(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.pdf",width = 3.5, height = 6.4,
#          antialias = "none", family = "Arial",  fallback_resolution = 600)


#pdf(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.pdf",width = 3.5, height = 6.2)
pheatmap::pheatmap(
  
  cluster_rlog_Sp[ , rownames(df1)],
  
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = F,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = F,
  border_color = "black",
  display_numbers = F,
  #clustering_method = "average", 
  #clustering_distance_rows = "maximum",
  cutree_cols = 2,
  #cutree_rows = 2,
  gaps_row = c( length(intersect(BregDown$Genes, SpleenDownNP$gene))-2,
                length(intersect(BregDown$Genes, SpleenDownNP$gene)) + length(intersect(BregUp$Genes, SpleenUpNP$gene)) -2,
                length(cluster_rlog_Sp$gene) -2
               ),
  gaps_col = 6,
  treeheight_col = 25,
  
  annotation_col = df1["Tx"],
  annotation_row = cluster_rlog_Sp["Expression in Bregs"],
  annotation_colors = anno_colors, 
  
  fontsize = 6.5,
  scale = "row", 
  fontsize_row = 5,
  cellwidth = 5,
  cellheight = 5,
  main = "Differentially Expressed\nBreg Markers"
)

dev.off()


#Spleen Volcano plot ----
Sp_deg_data = Spleen_NPvPBS3

keyvals <- ifelse(
  (Sp_deg_data$log2FoldChange < -LFc  & Sp_deg_data$padj < pc ), 'grey',
  ifelse((Sp_deg_data$log2FoldChange > LFc & Sp_deg_data$padj < pc), '#89CFF0',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'black'] <- 'Not DEG'
names(keyvals)[keyvals == '#89CFF0'] <- 'OVA NP'
names(keyvals)[keyvals == 'grey'] <- 'PBS'


TopDEGsSp = Sp_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, abs(log2FoldChange) > LFc) %>%
  arrange(desc(abs(log2FoldChange)))

TopDEGsUpSp = Sp_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, log2FoldChange > LFc) %>%
  arrange(desc(abs(log2FoldChange)))

TopDEGsDownSp = Sp_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, log2FoldChange < -LFc) %>%
  arrange(desc(abs(log2FoldChange)))

head(TopDEGsDownSp)

c(rownames(TopDEGsUpSp)[1:10],rownames(TopDEGsDownSp)[1:10])

#png(file = "SpleenBcells_VolcanoPlot_top20genes.png",units = "in",width = 5, height = 4.5, res = 400)
#pdf(file = "SpleenBcells_VolcanoPlot_top20genes.pdf",width = 6, height = 3.5)
EnhancedVolcano(Sp_deg_data,
                lab = rownames(Sp_deg_data),
                x = 'log2FoldChange',
                y = 'padj',
                axisLabSize = 7,
                titleLabSize = 5,
                subtitleLabSize = 5,
                captionLabSize = 1,
                pointSize = 1,
                labSize = 1.7,
                legendLabSize = 8,
                legendIconSize = 4,
                borderWidth = 0,
                vlineWidth = 0.2,
                hlineWidth = 0.2,
                gridlines.major = F,
                gridlines.minor = F,
                
                
                pCutoff = pc,
                FCcutoff = LFc,
                
                colCustom = keyvals,
                #selectLab = rownames(TopDEGs)[1:25],
                selectLab = c(rownames(TopDEGsUpSp)[1:10],rownames(TopDEGsDownSp)[1:10]),
                boxedLabels = T,
                drawConnectors = T,
                max.overlaps = 100,
                typeConnectors = "open", #(“open”, “closed”), ‘endsConnectors’ (“last”, “first”, “both”), 
                endsConnectors = "first",
                legendPosition = "top",
                title = NULL,
                caption = paste("Total DEGs", length(TopDEGsSp$log2FoldChange)),
                subtitle = "Spleen") + theme_prism(base_size = 12,base_fontface = "bold")+ 
  theme(axis.ticks = element_line(linewidth = 0.5), legend.position = "top", title = element_text(),
        legend.text = element_text(size = 8)
  )  + xlim(-1.3,2.4) #+ ylim(-1, 5)
#ggsave(file = "AllergyScaf_PBSLHCh4_VolcanoPlot.png", units = "in",width = 3, height = 3, res = 400)
dev.off()

#PP Volcano plot ----
PP_deg_data = PP_NPvPBS3

keyvals1 <- ifelse(
  (PP_deg_data$log2FoldChange < -LFc  & PP_deg_data$padj < pc ), 'grey',
  ifelse((PP_deg_data$log2FoldChange > LFc & PP_deg_data$padj < pc), '#89CFF0',
         'gray30'))
keyvals1[is.na(keyvals1)] <- 'gray30'
names(keyvals1)[keyvals1 == 'gray30'] <- 'Not DEG'
names(keyvals1)[keyvals1 == '#89CFF0'] <- 'OVA NP'
names(keyvals1)[keyvals1 == 'grey'] <- 'PBS'


TopDEGsPP=PP_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, abs(log2FoldChange) > LFc) %>%
  arrange(desc(abs(log2FoldChange)))

TopDEGsUpPP=PP_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, log2FoldChange > LFc) %>%
  arrange(desc(abs(log2FoldChange)))

TopDEGsDownPP=PP_deg_data %>%
  as.data.frame %>%
  filter(padj < pc, log2FoldChange < -LFc) %>%
  arrange(desc(abs(log2FoldChange)))

head(TopDEGsDownPP)

c(rownames(TopDEGsUpPP)[1:10],rownames(TopDEGsDownPP)[1:10])

#png(file = "PPBcells_VolcanoPlot_BregMarkers.png",units = "in",width = 5, height = 4.5, res = 400)

#pdf(file = "PPBcells_VolcanoPlot_BregMarkers.pdf",width = 4, height = 3.5)
EnhancedVolcano(PP_deg_data,
                lab = rownames(PP_deg_data),
                x = 'log2FoldChange',
                y = 'padj',
                axisLabSize = 7,
                titleLabSize = 5,
                subtitleLabSize = 5,
                captionLabSize = 2,
                pointSize = 1.5,
                labSize = 2,
                legendLabSize = 10,
                legendIconSize = 4,
                borderWidth = 0,
                vlineWidth = 0.2,
                hlineWidth = 0.2,
                gridlines.major = F,
                gridlines.minor = F,
                
                
                pCutoff = pc,
                FCcutoff = LFc,
                
                colCustom = keyvals1,
                #selectLab = rownames(TopDEGsPP)[1:25],
                #selectLab = c(rownames(TopDEGsUpPP)[1:10],rownames(TopDEGsDownPP)[1:10]),
                selectLab = c(
                              intersect(BregUp$Genes, PPUpNP$gene),
                              intersect(BregUp$Genes, PPDownNP$gene),
                              intersect(BregDown$Genes, PPDownNP$gene),
                              intersect(BregDown$Genes, PPUpNP$gene)
                              ),
                boxedLabels = T,
                drawConnectors = T,
                max.overlaps = 100,
                typeConnectors = "open", #(“open”, “closed”), ‘endsConnectors’ (“last”, “first”, “both”), 
                endsConnectors = "first",
                legendPosition = "top",
                title = NULL,
                caption = paste("Total DEGs", length(TopDEGsPP$log2FoldChange)),
                subtitle = "Peyer's Patches") + theme_prism(base_size = 12,base_fontface = "bold")+ 
  theme(axis.ticks = element_line(linewidth = 0.5), legend.position = "top", title = element_text(),
        legend.text = element_text(size = 9)
  )  + xlim(-1,2) + ylim(-5, 39)
#ggsave(file = "AllergyScaf_PBSLHCh4_VolcanoPlot.png", units = "in",width = 3, height = 3, res = 400)
dev.off()


# Functional Enrichment ----



#Spleen NP vPBS ----
DE_Sp_PBSvNP_EA = as.data.frame(Spleen_NPvPBS3)
DE_Sp_PBSvNP_EA$symbol = rownames(DE_Sp_PBSvNP_EA)
head(DE_Sp_PBSvNP_EA)

head(gene_info)
res_ids_sp <- inner_join(DE_Sp_PBSvNP_EA, gene_info, by=c("symbol"="external_gene_name")) 
head(res_ids_sp)

geneListDEGsSp = list(
  
  "NP" = res_ids_sp[which(res_ids_sp$log2FoldChange > 1  & res_ids_sp$padj < pc),"gene_id"],
  
  "PBS" = res_ids_sp[which(res_ids_sp$log2FoldChange < -0.25 & res_ids_sp$padj < pc),"gene_id"]
  
)


ck <- compareCluster(geneCluster = geneListDEGsSp, fun = enrichGO, 
                     universe      = data$gene_id,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE
)
#ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENSEMBL")
head(ck)

view(ck@compareClusterResult)

d <- godata('org.Mm.eg.db', ont="BP")
ck <- pairwise_termsim(ck, method="JC", semData = d,
                       showCategory = 500)



trace("treeplot.compareClusterResult", edit=TRUE, where = enrichplot::treeplot)

# edits: 
# line 142: dotdata$GeneRatio <- parse_ratio(dotdata$GeneRatio)
# line 146 size = "GeneRatio
# line 147: pwidth = 0.1
# remove + coord_equal() from last line

#png(file = "GOenrich_Spleen.png",units = "in",width = 14, height = 6, res = 400)
#cairo_pdf(file = "GOenrich_Spleen_dotplot.pdf",width = 14, height = 6)

treeplot(ck,
        #filter(ck, Count > 20),
         showCategory = 15, 
         fontsize = 4, cex_category = 0.9,
         #label_format_tiplab = 55,
         cluster.params = list(n = 6, 
                               label_words_n = 6, 
                               label_format = 35,
                               method = "mcquitty"
                               #color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
         ),
         offset.params = list(bar_tree = rel(1), tiplab = rel(1), extend = 0.2, hexpand = 0.15),
         hilight.params = list(hilight = T, align = "left"), 
        clusterPanel.params = list(clusterPanel = "dotplot", legend_n = 3, size = "GeneRatio")) + ggtitle("Spleen") #+theme_prism()
dev.off()


#PP NP v PBS----
DE_PP_PBSvNP_EA = as.data.frame(PP_NPvPBS3)
DE_PP_PBSvNP_EA$symbol = rownames(DE_PP_PBSvNP_EA)

res_ids_pp <- inner_join(DE_PP_PBSvNP_EA, gene_info, by=c("symbol"="external_gene_name")) 

geneListDEGsPP = list(
  
  "NP" = res_ids_pp[which(res_ids_pp$log2FoldChange > 0  & res_ids_pp$padj < pc),"gene_id"],
  
  "PBS" = res_ids_pp[which(res_ids_pp$log2FoldChange < 0  & res_ids_pp$padj < pc),"gene_id"]
  
)


ck4 <- compareCluster(geneCluster  = geneListDEGsPP, fun = enrichGO, 
                     universe      = data$gene_id,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "fdr",
                     pvalueCutoff  = 0.1,
                     qvalueCutoff  = 0.1,
                     readable      = TRUE
)
head(ck4) 

#d <- godata('org.Mm.eg.db', ont="BP")
ck4 <- pairwise_termsim(ck4, method="JC", semData = d,
                       showCategory = 500)
#cairo_pdf(file = "GOenrich_PP_dotplot.pdf",width = 14, height = 6)
treeplot(ck4,
         #filter(ck, Count > 20),
         showCategory = 15, 
         fontsize = 4, cex_category = 0.9,
         #label_format_tiplab = 55,
         cluster.params = list(n = 6, 
                               label_words_n = 6, 
                               label_format = 35,
                               method = "mcquitty"
                               #color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
         ),
         offset.params = list(bar_tree = rel(1), tiplab = rel(1), extend = 0.2, hexpand = 0.15),
         hilight.params = list(hilight = T, align = "left"), 
         clusterPanel.params = list(clusterPanel = "dotplot", legend_n = 3, size = "GeneRatio")) + ggtitle("PP")
dev.off()



#plot individual genes ----

#plotCounts(dds_full_Spleen, gene = "Fcer2a", intgroup = c("Tx"))

#SI FIG 11 ----
#plot CD23 expression ----
#Spleen CD23 ----
fiss <- plotCounts(dds_full_Spleen, 
                   "Fcer2a", normalized = T,
                   intgroup = c("Tx","Tissue","Sample Label"), returnData = T)

fiss = fiss[fiss$Tx %in% c("PBS", "OVA NPs"),]
fiss$Tx = factor(fiss$Tx, levels = c("PBS", "OVA NPs") )
fiss
Spleen_filt1["Fcer2a",]

Spleen_NPvPBS3["Fcer2a",]
Spleen_sOVAvPBS3["Fcer2a",]
Spleen_NPvsOVA3["Fcer2a",]



stat.test <- tibble::tribble(
  ~group1, ~group2,   ~p.adj,
  "PBS",     "OVA NPs", round(Spleen_NPvPBS3["Fcer2a","padj"],5)
)
stat.test = stat.test %>%
  add_x_position()%>%
  add_significance()

stat.test

Spleen_CD23_p = ggplot(fiss,
       aes(x = Tx, y = log2(count))) + 
  geom_violin(aes(x = Tx), scale = "area", trim = F ) +
  stat_summary( fun=mean, geom= "crossbar", color = "red", width = 0.1, size = 0.3) +
  stat_summary(
    fun.min = function(x) mean(x) - sd(x), 
    fun.max = function(x) mean(x) + sd(x), 
    geom = "errorbar", width = 0.1,
    color = "red", size = 0.3
  ) +
  geom_point(size = 0.5, alpha = 0.5) +
  #scale_y_log10() + 
  #scale_y_continuous(trans = "log2", breaks = c(0, 2, 2^2, 2^3, 2^4, 2^5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)) , breaks = seq(0, 20, 1), limits = c(11.7,15.5))+
  #scale_color_manual(values = c("purple3", "orange1")) +
  labs(x = "Spleen B Cells", title = "Pre-OFC 1", subtitle = "Fcer2a (CD23)") +
  theme_prism(base_family = "sans",  base_size = 10)  +
  stat_pvalue_manual(
    stat.test[1,],  y.position =14.8, 
    label = "{p.adj.signif}\n{p.adj}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) 
  
Spleen_CD23_p

#PP CD23----
fiss1 <- plotCounts(dds_NP_PP, 
                   "Fcer2a", normalized = T,
                   intgroup = c("Tx","Tissue","Sample Label"), returnData = T)


#fiss1 = PP_info[c("Tx","Tissue","Sample Label")]
#fiss1$count = as.numeric(PP["Fcer2a",])
class(fiss1$count)
fiss1
fiss1$Tx

PP_NPvPBS3["Fcer2a",]


stat.test1 <- tibble::tribble(
  ~group1, ~group2,   ~p.adj, ~p.adj.label,
  "PBS",     "OVA NPs", PP_NPvPBS3["Fcer2a","padj"], scientific(PP_NPvPBS3["Fcer2a","padj"]),
)
stat.test1 = stat.test1 %>%
  add_x_position()%>%
  add_significance()

stat.test1

PP_CD23_p = ggplot(fiss1,
       aes(x = Tx, y = log2(count))) + 
  geom_violin(aes(x = Tx), trim = F) +
  stat_summary( fun=mean, geom= "crossbar", color = "red", width = 0.1, size = 0.3) +
  stat_summary(
    fun.min = function(x) mean(x) - sd(x), 
    fun.max = function(x) mean(x) + sd(x), 
    geom = "errorbar", width = 0.1,
    color = "red",  size = 0.3
  ) +
  geom_point(size = 0.5,alpha = 0.5) +
  #scale_y_log10() + 
  #scale_y_continuous(trans = "log2", breaks = c(0, 2, 2^2, 2^3, 2^4, 2^5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)) ,  breaks = seq(0, 20, 1), limits = c(11.7,15.5))+
  #scale_color_manual(values = c("purple3", "orange1")) +
  labs(x = "PP B Cells", title = "Pre-OFC 1", subtitle = "Fcer2a (CD23)") +
  theme_prism(base_family = "sans", base_size = 10)  +
  stat_pvalue_manual(
    stat.test1[1,],  y.position =14.9, 
    label = "{p.adj.signif}\n{p.adj.label}",tip.length = 0.01, bracket.size = 0.3, #bracket.nudge.y = 700000,
    lineheight = 0.1, size = 2
  ) 

PP_CD23_p

pdf("Fcer2a_Bcells_SpleenPP_vlnplot.pdf", width = 4.4, 2.8)
ggarrange(S