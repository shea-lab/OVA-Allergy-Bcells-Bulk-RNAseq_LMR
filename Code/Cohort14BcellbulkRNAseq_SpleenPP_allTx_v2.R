

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
library(ggkegg)
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

#194 Breg Genes
length(unique(c(BregUp$Genes, BregDown$Genes)))

#42 OVA NP DEGs are Breg Markers
length(
unique(
c(
intersect(BregUp$Genes, SpleenUpNP$gene),
intersect(BregDown$Genes, SpleenDownNP$gene)
)))

39/186

length(unique(intersect(BregUp$Genes, SpleenDownNP$gene)))



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
  intersect(BregDown$Genes, SpleenDownNP$gene),
  intersect(BregDown$Genes, SpleenUpNP$gene)
)

n_occur <- data.frame(table(BregGenesMap))
n_occur[n_occur$Freq > 1,]

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
length(BregMarkersDEGsSp)

length(BregGenesMap)

# Remove duplicate genes
BregGenesMap = unique(BregGenesMap)
length(BregGenesMap)

# Save information about Breg Markers that are DEGs in a seperate csv file
#write.csv( Bcellmarkers[which(Bcellmarkers$Genes %in% BregGenesMap & Bcellmarkers$Subset == "Regulatory"), ], 
#           "BregMarkersDEGsfromR5.csv")

# Dataframe with only the Breg markers that are DEGs
BregmarkersDEG = read.csv("BregMarkersDEGsfromR5.csv")
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


view(cluster_rlog_Sp[intersect(BregUp$Genes, SpleenUpNP$gene),])
view(cluster_rlog_Sp[intersect(BregDown$Genes, SpleenDownNP$gene),])

view(cluster_rlog_Sp[intersect(BregUp$Genes, SpleenDownNP$gene),])



# Heatmap of Breg Markers DEGs in Spleen ----

#png(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.png",units = "in",width = 3.5, height = 6.2, res = 400)
#cairo_pdf(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.pdf",width = 3.5, height = 6.4,
#          antialias = "none", family = "Arial",  fallback_resolution = 600)


pdf(file = "SpleenBcells_Heatmap_DEGs_BregMarkers.pdf",width = 3.5, height = 6.2)
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



#Spleen NP v PBS ----
DE_Sp_PBSvNP_EA = as.data.frame(Spleen_NPvPBS3)
DE_Sp_PBSvNP_EA$symbol = rownames(DE_Sp_PBSvNP_EA)
head(DE_Sp_PBSvNP_EA)

head(gene_info)
res_ids_sp <- inner_join(DE_Sp_PBSvNP_EA, gene_info, by=c("symbol"="external_gene_name")) 
head(res_ids_sp)

geneListDEGsSp = list(
  
  "NP" = res_ids_sp[which(res_ids_sp$log2FoldChange > 0  & res_ids_sp$padj < pc),"gene_id"],
  
  "PBS" = res_ids_sp[which(res_ids_sp$log2FoldChange < 0 & res_ids_sp$padj < pc),"gene_id"]
  
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
# packageVersion("ggnewscale") #'0.5.2' if not: 
# devtools::install_github("eliocamp/ggnewscale@v0.5.2")
# line 142: dotdata$GeneRatio <- parse_ratio(dotdata$GeneRatio)
# line 146 size = "GeneRatio
# line 147: pwidth = 0.1, offset = 0, text.size = 6
# remove + coord_equal() from last line

#png(file = "GOenrich_Spleen.png",units = "in",width = 14, height = 6, res = 400)
#cairo_pdf(file = "GOenrich_Spleen_dotplot2.pdf",width = 14, height = 12)

treeplot(ck,
        #filter(ck, Count > 20),
         showCategory = 15, 
         fontsize = 6, cex_category = 6,
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

 #gseKEGG Sp NPvPBS ----
head(res_ids_sp)
Sp_geneList =  res_ids_sp %>%
  #filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  dplyr::select(entrezgene_id, log2FoldChange)
head(Sp_geneList)

Sp_geneList = Sp_geneList[!Sp_geneList$"entrezgene_id" == ".", ]

Sp_geneList <- deframe(Sp_geneList)
length(Sp_geneList)
Sp_geneList = na.omit(Sp_geneList)
length(Sp_geneList)
view(Sp_geneList)
length(unique(names(Sp_geneList)))
length((names(Sp_geneList)))


Sp_kegg <- gseKEGG(geneList     = Sp_geneList,
               organism     = 'mmu',
               #minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(Sp_kegg)
view(Sp_kegg)

dotplot(Sp_kegg, color = "NES", size = "p.adjust",
        showCategory = 30)
getwd()
png(file = "Sp_KEGG_GSEA.png",units = "in",width = 10, height = 12, res = 400)
dotplot(Sp_kegg, x = "NES", color = "p.adjust", size = "GeneRatio" ,orderBy = "NES",
        showCategory = 30, label_format = 60, font.size = 18) + 
  labs( title = "KEGG GSEA") + 
  geom_vline(xintercept = 0, size = 1) +
  theme(title = element_text(size = 20 ),
        plot.title.position = "panel"
  ) +
  theme_prism(base_line_size = 1)
dev.off()

Sp_kegg@result[1, c("Description")]
pdf("p53 signaling pathway_KEGG_GSEA_DEGs.pdf", width = 7, height = 5)
rawValue(Sp_geneList, Sp_kegg@result$ID[1], auto_add=TRUE ) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red4",
                       limits = c(-2,2)) +
  labs(fill = "logFc")


rawMap(Sp_kegg, Sp_kegg@result$ID[1])
rawValue(Sp_geneList, Sp_kegg@result$ID[1], auto_add=TRUE)
length(unique(Sp_geneList))
length(Sp_geneList)

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
count.data["Il10",]
plotCounts(dds_full_Spleen, 
           "Il10", normalized = T,
           intgroup = c("Tx","Tissue","Sample Label"), returnData = T)

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
ggarrange(Spleen_CD23_p, PP_CD23_p,ncol = 2, labels = "AUTO")
dev.off()

getwd()

# Session Info ----
# devtools::session_info()
# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS 15.2
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Detroit
# date     2025-07-26
# rstudio  2023.03.1+446 Cherry Blossom (desktop)
# pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                    1.4-5      2016-07-21 [2] CRAN (R 4.3.0)
# annotate                 1.78.0     2023-04-25 [2] Bioconductor
# AnnotationDbi          * 1.62.0     2023-04-25 [2] Bioconductor
# AnnotationFilter       * 1.24.0     2023-04-25 [2] Bioconductor
# AnnotationHub          * 3.8.0      2023-04-25 [2] Bioconductor
# ape                      5.7-1      2023-03-13 [2] CRAN (R 4.3.0)
# apeglm                 * 1.21.0     2022-12-20 [2] Bioconductor
# aplot                    0.1.10     2023-03-08 [2] CRAN (R 4.3.0)
# babelgene                22.9       2022-09-29 [2] CRAN (R 4.3.0)
# backports                1.4.1      2021-12-13 [2] CRAN (R 4.3.0)
# base64enc                0.1-3      2015-07-28 [2] CRAN (R 4.3.0)
# bbmle                    1.0.25     2022-05-11 [2] CRAN (R 4.3.0)
# bdsmatrix                1.3-6      2022-06-03 [2] CRAN (R 4.3.0)
# Biobase                * 2.59.0     2022-12-20 [2] Bioconductor
# BiocFileCache          * 2.8.0      2023-04-25 [2] Bioconductor
# BiocGenerics           * 0.46.0     2023-04-25 [2] Bioconductor
# BiocIO                   1.10.0     2023-04-25 [2] Bioconductor
# BiocManager              1.30.20    2023-02-24 [2] CRAN (R 4.3.0)
# BiocNeighbors            1.17.1     2022-12-20 [2] Bioconductor
# BiocParallel             1.33.11    2023-03-24 [2] Bioconductor
# BiocVersion              3.17.1     2022-12-20 [2] Bioconductor
# biomaRt                * 2.56.0     2023-04-25 [2] Bioconductor
# Biostrings               2.67.2     2023-04-19 [2] Bioconductor
# bit                      4.0.5      2022-11-15 [2] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [2] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [2] CRAN (R 4.3.0)
# blob                     1.2.4      2023-03-17 [2] CRAN (R 4.3.0)
# boot                   * 1.3-28.1   2022-11-22 [2] CRAN (R 4.3.0)
# broom                    1.0.4      2023-03-11 [2] CRAN (R 4.3.0)
# cachem                   1.0.8      2023-05-01 [2] CRAN (R 4.3.0)
# callr                    3.7.3      2022-11-02 [2] CRAN (R 4.3.0)
# car                    * 3.1-2      2023-03-30 [2] CRAN (R 4.3.0)
# carData                * 3.0-5      2022-01-06 [2] CRAN (R 4.3.0)
# caret                  * 6.0-94     2023-03-21 [2] CRAN (R 4.3.0)
# caTools                  1.18.2     2021-03-28 [2] CRAN (R 4.3.0)
# CellChat                 1.6.1      2023-05-21 [2] Github (sqjin/CellChat@268b447)
# checkmate                2.2.0      2023-04-27 [2] CRAN (R 4.3.0)
# circlize                 0.4.15     2022-05-10 [2] CRAN (R 4.3.0)
# class                    7.3-21     2023-01-23 [2] CRAN (R 4.3.0)
# cli                      3.6.1      2023-03-23 [2] CRAN (R 4.3.0)
# clue                     0.3-64     2023-01-31 [2] CRAN (R 4.3.0)
# cluster                  2.1.4      2022-08-22 [2] CRAN (R 4.3.0)
# clusterProfiler        * 4.8.0      2023-04-25 [2] Bioconductor
# coda                     0.19-4     2020-09-30 [2] CRAN (R 4.3.0)
# codetools                0.2-19     2023-02-01 [2] CRAN (R 4.3.0)
# coefplot               * 1.2.8      2022-01-14 [2] CRAN (R 4.3.0)
# colorspace             * 2.1-0      2023-01-23 [2] CRAN (R 4.3.0)
# ComplexHeatmap           2.15.4     2023-05-21 [2] Github (jokergoo/ComplexHeatmap@ae0ec42)
# ConsensusClusterPlus     1.64.0     2023-05-08 [1] Bioconductor
# corpcor                  1.6.10     2021-09-16 [2] CRAN (R 4.3.0)
# cowplot                * 1.1.1      2020-12-30 [2] CRAN (R 4.3.0)
# crayon                   1.5.2      2022-09-29 [2] CRAN (R 4.3.0)
# curl                     5.0.0      2023-01-12 [2] CRAN (R 4.3.0)
# data.table               1.14.8     2023-02-17 [2] CRAN (R 4.3.0)
# datawizard               1.1.0      2025-05-09 [1] CRAN (R 4.3.3)
# DBI                      1.1.3      2022-06-18 [2] CRAN (R 4.3.0)
# dbplyr                 * 2.3.2      2023-03-21 [2] CRAN (R 4.3.0)
# DEGreport              * 1.36.0     2023-05-08 [1] Bioconductor
# DelayedArray             0.25.0     2022-12-20 [2] Bioconductor
# dendextend               1.17.1     2023-03-25 [2] CRAN (R 4.3.0)
# DEoptimR                 1.0-13     2023-05-02 [2] CRAN (R 4.3.0)
# DESeq2                 * 1.40.2     2023-07-02 [1] Bioconductor
# devtools               * 2.4.5      2022-10-11 [2] CRAN (R 4.3.0)
# digest                   0.6.31     2022-12-11 [2] CRAN (R 4.3.0)
# doParallel               1.0.17     2022-02-07 [2] CRAN (R 4.3.0)
# DOSE                   * 3.26.0     2023-04-25 [2] Bioconductor
# downloader               0.4        2015-07-09 [2] CRAN (R 4.3.0)
# dplyr                  * 1.1.2      2023-04-20 [2] CRAN (R 4.3.0)
# edgeR                  * 3.42.2     2023-05-08 [2] Bioconductor
# ellipse                  0.4.5      2023-04-05 [2] CRAN (R 4.3.0)
# ellipsis                 0.3.2      2021-04-29 [2] CRAN (R 4.3.0)
# emdbook                  1.3.12     2020-02-19 [2] CRAN (R 4.3.0)
# EnhancedVolcano        * 1.18.0     2023-05-08 [1] Bioconductor
# enrichplot             * 1.20.0     2023-04-25 [2] Bioconductor
# enrichR                * 3.2        2023-04-14 [2] CRAN (R 4.3.0)
# ensembldb              * 2.23.2     2023-02-15 [2] Bioconductor
# evaluate                 0.20       2023-01-17 [2] CRAN (R 4.3.0)
# factoextra             * 1.0.7      2020-04-01 [2] CRAN (R 4.3.0)
# fansi                    1.0.4      2023-01-22 [2] CRAN (R 4.3.0)
# farver                   2.1.1      2022-07-06 [2] CRAN (R 4.3.0)
# fastmap                  1.1.1      2023-02-24 [2] CRAN (R 4.3.0)
# fastmatch                1.1-3      2021-07-23 [2] CRAN (R 4.3.0)
# fgsea                    1.25.0     2022-12-20 [2] Bioconductor
# filelock                 1.0.2      2018-10-05 [2] CRAN (R 4.3.0)
# fmsb                   * 0.7.5      2023-01-05 [2] CRAN (R 4.3.0)
# FNN                      1.1.3.2    2023-03-20 [2] CRAN (R 4.3.0)
# forcats                * 1.0.0      2023-01-29 [2] CRAN (R 4.3.0)
# foreach                  1.5.2      2022-02-02 [2] CRAN (R 4.3.0)
# foreign                  0.8-84     2022-12-06 [2] CRAN (R 4.3.0)
# Formula                  1.2-5      2023-02-24 [2] CRAN (R 4.3.0)
# fs                       1.6.4      2024-04-25 [2] CRAN (R 4.3.1)
# future                   1.32.0     2023-03-07 [2] CRAN (R 4.3.0)
# future.apply             1.10.0     2022-11-05 [2] CRAN (R 4.3.0)
# genefilter             * 1.82.1     2023-05-08 [2] Bioconductor
# generics                 0.1.3      2022-07-05 [2] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.36.0     2023-04-25 [2] Bioconductor
# GenomeInfoDbData         1.2.10     2023-04-24 [2] Bioconductor
# GenomicAlignments        1.35.1     2023-03-24 [2] Bioconductor
# GenomicFeatures        * 1.51.4     2022-12-23 [2] Bioconductor
# GenomicRanges          * 1.51.4     2022-12-20 [2] Bioconductor
# GetoptLong               1.0.5      2020-12-15 [2] CRAN (R 4.3.0)
# ggalluvial               0.12.5     2023-02-22 [2] CRAN (R 4.3.0)
# ggbiplot               * 0.55       2023-04-24 [2] Github (vqv/ggbiplot@7325e88)
# ggdendro                 0.1.23     2022-02-16 [2] CRAN (R 4.3.0)
# ggforce                * 0.4.1      2022-10-04 [2] CRAN (R 4.3.0)
# ggfun                    0.0.9      2022-11-21 [2] CRAN (R 4.3.0)
# ggnetwork                0.5.12     2023-03-06 [2] CRAN (R 4.3.0)
# ggnewscale               0.4.8      2022-10-06 [2] CRAN (R 4.3.0)
# ggplot2                * 3.4.2      2023-04-03 [2] CRAN (R 4.3.0)
# ggplotify                0.1.0      2021-09-02 [2] CRAN (R 4.3.0)
# ggprism                * 1.0.4      2022-11-04 [2] CRAN (R 4.3.0)
# ggpubr                 * 0.6.0      2023-02-10 [2] CRAN (R 4.3.0)
# ggraph                   2.1.0      2022-10-09 [2] CRAN (R 4.3.0)
# ggrepel                * 0.9.3      2023-02-03 [2] CRAN (R 4.3.0)
# ggsignif                 0.6.4      2022-10-13 [2] CRAN (R 4.3.0)
# ggtree                   3.8.0      2023-04-25 [2] Bioconductor
# ggtreeExtra            * 1.10.0     2023-05-08 [1] Bioconductor
# ggvenn                 * 0.1.10     2023-03-31 [2] CRAN (R 4.3.0)
# ggwordcloud              0.5.0      2019-06-02 [2] CRAN (R 4.3.0)
# glmnet                 * 4.1-7      2023-03-23 [2] CRAN (R 4.3.0)
# GlobalOptions            0.1.2      2020-06-10 [2] CRAN (R 4.3.0)
# globals                  0.16.2     2022-11-21 [2] CRAN (R 4.3.0)
# glue                     1.6.2      2022-02-24 [2] CRAN (R 4.3.0)
# GO.db                    3.17.0     2023-04-24 [2] Bioconductor
# GOSemSim               * 2.25.0     2022-12-20 [2] Bioconductor
# gower                    1.0.1      2022-12-22 [2] CRAN (R 4.3.0)
# gplots                 * 3.1.3      2022-04-25 [2] CRAN (R 4.3.0)
# graph                    1.77.2     2023-03-24 [2] Bioconductor
# graphlayouts             1.0.0      2023-05-01 [2] CRAN (R 4.3.0)
# gridBase                 0.4-7      2014-02-24 [2] CRAN (R 4.3.0)
# gridExtra              * 2.3        2017-09-09 [2] CRAN (R 4.3.0)
# gridGraphics             0.5-1      2020-12-13 [2] CRAN (R 4.3.0)
# GSEAmining             * 1.10.0     2023-04-25 [2] Bioconductor
# gson                     0.1.0      2023-03-07 [2] CRAN (R 4.3.0)
# gtable                   0.3.3      2023-03-21 [2] CRAN (R 4.3.0)
# gtools                   3.9.4      2022-11-27 [2] CRAN (R 4.3.0)
# hardhat                  1.3.0      2023-03-30 [2] CRAN (R 4.3.0)
# HDO.db                   0.99.1     2023-04-24 [2] Bioconductor
# Hmisc                  * 5.0-1      2023-03-08 [2] CRAN (R 4.3.0)
# hms                      1.1.3      2023-03-21 [2] CRAN (R 4.3.0)
# htmlTable                2.4.1      2022-07-07 [2] CRAN (R 4.3.0)
# htmltools                0.5.5      2023-03-23 [2] CRAN (R 4.3.0)
# htmlwidgets              1.6.2      2023-03-17 [2] CRAN (R 4.3.0)
# httpuv                   1.6.9      2023-02-14 [2] CRAN (R 4.3.0)
# httr                     1.4.5      2023-02-24 [2] CRAN (R 4.3.0)
# igraph                   1.4.2      2023-04-07 [2] CRAN (R 4.3.0)
# insight                  1.2.0      2025-04-22 [1] CRAN (R 4.3.3)
# interactiveDisplayBase   1.38.0     2023-04-25 [2] Bioconductor
# ipred                    0.9-14     2023-03-09 [2] CRAN (R 4.3.0)
# IRanges                * 2.33.1     2023-04-19 [2] Bioconductor
# irlba                    2.3.5.1    2022-10-03 [2] CRAN (R 4.3.0)
# iterators                1.0.14     2022-02-05 [2] CRAN (R 4.3.0)
# janeaustenr              1.0.0      2022-08-26 [2] CRAN (R 4.3.0)
# jsonlite                 1.8.4      2022-12-06 [2] CRAN (R 4.3.0)
# KEGGgraph                1.60.0     2023-04-25 [2] Bioconductor
# KEGGREST                 1.40.0     2023-04-25 [2] Bioconductor
# KernSmooth               2.23-20    2021-05-03 [2] CRAN (R 4.3.0)
# knitr                    1.42       2023-01-25 [2] CRAN (R 4.3.0)
# later                    1.3.1      2023-05-02 [2] CRAN (R 4.3.0)
# lattice                * 0.21-8     2023-04-05 [2] CRAN (R 4.3.0)
# lava                     1.7.2.1    2023-02-27 [2] CRAN (R 4.3.0)
# lazyeval                 0.2.2      2019-03-15 [2] CRAN (R 4.3.0)
# lifecycle                1.0.3      2022-10-07 [2] CRAN (R 4.3.0)
# limma                  * 3.55.7     2023-04-19 [2] Bioconductor
# listenv                  0.9.0      2022-12-16 [2] CRAN (R 4.3.0)
# locfit                   1.5-9.7    2023-01-02 [2] CRAN (R 4.3.0)
# logging                  0.10-108   2019-07-14 [1] CRAN (R 4.3.3)
# lubridate              * 1.9.2      2023-02-10 [2] CRAN (R 4.3.0)
# magrittr               * 2.0.3      2022-03-30 [2] CRAN (R 4.3.0)
# MASS                   * 7.3-59     2023-04-21 [2] CRAN (R 4.3.0)
# Matrix                 * 1.5-4      2023-04-04 [2] CRAN (R 4.3.0)
# MatrixGenerics         * 1.12.0     2023-04-25 [2] Bioconductor
# matrixStats            * 1.3.0      2024-04-11 [1] CRAN (R 4.3.1)
# memoise                  2.0.1      2021-11-26 [2] CRAN (R 4.3.0)
# MeSHDbi                  1.36.0     2023-04-25 [2] Bioconductor
# mime                     0.12       2021-09-28 [2] CRAN (R 4.3.0)
# miniUI                   0.1.1.1    2018-05-18 [2] CRAN (R 4.3.0)
# minpack.lm             * 1.2-3      2023-01-26 [2] CRAN (R 4.3.0)
# mixOmics               * 6.24.0     2023-04-25 [2] Bioconductor (R 4.3.0)
# mnormt                   2.1.1      2022-09-26 [2] CRAN (R 4.3.0)
# ModelMetrics             1.2.2.2    2020-03-17 [2] CRAN (R 4.3.0)
# msigdbr                * 7.5.1      2022-03-30 [2] CRAN (R 4.3.0)
# munsell                  0.5.0      2018-06-12 [2] CRAN (R 4.3.0)
# mvtnorm                  1.1-3      2021-10-08 [2] CRAN (R 4.3.0)
# network                  1.18.1     2023-01-24 [2] CRAN (R 4.3.0)
# nlme                     3.1-162    2023-01-31 [2] CRAN (R 4.3.0)
# NMF                      0.26       2023-03-20 [2] CRAN (R 4.3.0)
# nnet                     7.3-18     2022-09-28 [2] CRAN (R 4.3.0)
# NormqPCR               * 1.46.0     2023-05-08 [2] Bioconductor
# numDeriv                 2016.8-1.1 2019-06-06 [2] CRAN (R 4.3.0)
# org.Hs.eg.db             3.17.0     2023-04-24 [2] Bioconductor
# org.Mm.eg.db           * 3.17.0     2023-04-24 [2] Bioconductor
# pacman                   0.5.1      2019-03-11 [2] CRAN (R 4.3.0)
# parallelly               1.35.0     2023-03-23 [2] CRAN (R 4.3.0)
# patchwork                1.1.2      2022-08-19 [2] CRAN (R 4.3.0)
# pathview               * 1.40.0     2023-04-25 [2] Bioconductor
# pbapply                  1.7-0      2023-01-13 [2] CRAN (R 4.3.0)
# performance              0.13.0     2025-01-15 [1] CRAN (R 4.3.3)
# pheatmap               * 1.0.12     2019-01-04 [2] CRAN (R 4.3.0)
# pillar                   1.9.0      2023-03-22 [2] CRAN (R 4.3.0)
# pkgbuild                 1.4.0      2022-11-27 [2] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [2] CRAN (R 4.3.0)
# pkgload                  1.3.2      2022-11-16 [2] CRAN (R 4.3.0)
# plyr                   * 1.8.8      2022-11-11 [2] CRAN (R 4.3.0)
# png                      0.1-8      2022-11-29 [2] CRAN (R 4.3.0)
# polyclip                 1.10-4     2022-10-20 [2] CRAN (R 4.3.0)
# prettyunits              1.1.1      2020-01-24 [2] CRAN (R 4.3.0)
# pROC                     1.18.0     2021-09-03 [2] CRAN (R 4.3.0)
# processx                 3.8.1      2023-04-18 [2] CRAN (R 4.3.0)
# prodlim                  2023.03.31 2023-04-02 [2] CRAN (R 4.3.0)
# profvis                  0.3.8      2023-05-02 [2] CRAN (R 4.3.0)
# progress                 1.2.2      2019-05-16 [2] CRAN (R 4.3.0)
# promises                 1.2.0.1    2021-02-11 [2] CRAN (R 4.3.0)
# ProtGenerics             1.32.0     2023-04-25 [2] Bioconductor
# ps                       1.7.5      2023-04-18 [2] CRAN (R 4.3.0)
# psych                    2.3.3      2023-03-18 [2] CRAN (R 4.3.0)
# purrr                  * 1.0.1      2023-01-10 [2] CRAN (R 4.3.0)
# qpcR                   * 1.4-1      2018-06-14 [1] CRAN (R 4.3.0)
# qvalue                   2.32.0     2023-04-25 [2] Bioconductor
# R6                       2.5.1      2021-08-19 [2] CRAN (R 4.3.0)
# randomForest           * 4.7-1.1    2022-05-23 [2] CRAN (R 4.3.0)
# RANN                     2.6.1      2019-01-08 [2] CRAN (R 4.3.0)
# rappdirs                 0.3.3      2021-01-31 [2] CRAN (R 4.3.0)
# rARPACK                  0.11-0     2016-03-10 [2] CRAN (R 4.3.0)
# RColorBrewer           * 1.1-3      2022-04-03 [2] CRAN (R 4.3.0)
# Rcpp                     1.0.10     2023-01-22 [2] CRAN (R 4.3.0)
# RcppAnnoy                0.0.20     2022-10-27 [2] CRAN (R 4.3.0)
# RcppHNSW                 0.4.1      2022-07-18 [2] CRAN (R 4.3.0)
# RCurl                    1.98-1.12  2023-03-27 [2] CRAN (R 4.3.0)
# rdist                  * 0.0.5      2020-05-04 [2] CRAN (R 4.3.0)
# ReadqPCR               * 1.46.0     2023-04-25 [2] Bioconductor
# readr                  * 2.1.4      2023-02-10 [2] CRAN (R 4.3.0)
# recipes                  1.0.6      2023-04-25 [2] CRAN (R 4.3.0)
# registry                 0.5-1      2019-03-05 [2] CRAN (R 4.3.0)
# remotes                  2.4.2      2021-11-30 [2] CRAN (R 4.3.0)
# reshape                  0.8.9      2022-04-12 [2] CRAN (R 4.3.0)
# reshape2               * 1.4.4      2020-04-09 [2] CRAN (R 4.3.0)
# restfulr                 0.0.15     2022-06-16 [2] CRAN (R 4.3.0)
# reticulate               1.28       2023-01-27 [2] CRAN (R 4.3.0)
# rgl                    * 1.3.1      2024-03-05 [1] CRAN (R 4.3.1)
# Rgraphviz                2.43.0     2022-12-20 [2] Bioconductor
# rjson                    0.2.21     2022-01-09 [2] CRAN (R 4.3.0)
# rlang                    1.1.1      2023-04-28 [2] CRAN (R 4.3.0)
# rmarkdown                2.21       2023-03-26 [2] CRAN (R 4.3.0)
# rngtools                 1.5.2      2021-09-20 [2] CRAN (R 4.3.0)
# robustbase             * 0.95-1     2023-03-29 [2] CRAN (R 4.3.0)
# ROCR                   * 1.0-11     2020-05-02 [2] CRAN (R 4.3.0)
# rpart                    4.1.19     2022-10-21 [2] CRAN (R 4.3.0)
# Rsamtools                2.15.2     2023-03-24 [2] Bioconductor
# RSpectra                 0.16-1     2022-04-24 [2] CRAN (R 4.3.0)
# RSQLite                  2.3.1      2023-04-03 [2] CRAN (R 4.3.0)
# rstatix                  0.7.2      2023-02-01 [2] CRAN (R 4.3.0)
# rstudioapi               0.14       2022-08-22 [2] CRAN (R 4.3.0)
# rtracklayer              1.59.1     2022-12-30 [2] Bioconductor
# Rtsne                    0.16       2022-04-17 [2] CRAN (R 4.3.0)
# S4Arrays                 1.2.1      2024-03-05 [1] Bioconductor 3.18 (R 4.3.2)
# S4Vectors              * 0.37.4     2023-03-24 [2] Bioconductor
# scales                 * 1.2.1      2022-08-20 [2] CRAN (R 4.3.0)
# scatterpie               0.1.9      2023-04-22 [2] CRAN (R 4.3.0)
# sessioninfo              1.2.2      2021-12-06 [2] CRAN (R 4.3.0)
# shadowtext               0.1.2      2022-04-22 [2] CRAN (R 4.3.0)
# shape                    1.4.6      2021-05-19 [2] CRAN (R 4.3.0)
# shiny                    1.7.4      2022-12-15 [2] CRAN (R 4.3.0)
# sjstats                  0.19.0     2024-05-14 [1] CRAN (R 4.3.3)
# sna                      2.7-1      2023-01-24 [2] CRAN (R 4.3.0)
# SnowballC                0.7.1      2023-04-25 [2] CRAN (R 4.3.0)
# statnet.common           4.8.0      2023-01-24 [2] CRAN (R 4.3.0)
# stringi                  1.7.12     2023-01-11 [2] CRAN (R 4.3.0)
# stringr                * 1.5.0      2022-12-02 [2] CRAN (R 4.3.0)
# SummarizedExperiment   * 1.29.1     2022-12-20 [2] Bioconductor
# survival                 3.5-5      2023-03-12 [2] CRAN (R 4.3.0)
# svglite                  2.1.3      2023-12-08 [2] CRAN (R 4.3.1)
# systemfonts              1.0.4      2022-02-11 [2] CRAN (R 4.3.0)
# tibble                 * 3.2.1      2023-03-20 [2] CRAN (R 4.3.0)
# tidygraph                1.2.3      2023-02-01 [2] CRAN (R 4.3.0)
# tidyr                  * 1.3.0      2023-01-24 [2] CRAN (R 4.3.0)
# tidyselect               1.2.0      2022-10-10 [2] CRAN (R 4.3.0)
# tidytext                 0.4.1      2023-01-07 [2] CRAN (R 4.3.0)
# tidytree                 0.4.2      2022-12-18 [2] CRAN (R 4.3.0)
# tidyverse              * 2.0.0      2023-02-22 [2] CRAN (R 4.3.0)
# timechange               0.2.0      2023-01-11 [2] CRAN (R 4.3.0)
# timeDate                 4022.108   2023-01-07 [2] CRAN (R 4.3.0)
# tokenizers               0.3.0      2022-12-22 [2] CRAN (R 4.3.0)
# treeio                   1.24.0     2023-04-25 [2] Bioconductor
# tweenr                   2.0.2      2022-09-06 [2] CRAN (R 4.3.0)
# tzdb                     0.3.0      2022-03-28 [2] CRAN (R 4.3.0)
# urlchecker               1.0.1      2021-11-30 [2] CRAN (R 4.3.0)
# useful                   1.2.6      2018-10-08 [2] CRAN (R 4.3.0)
# usethis                * 2.1.6      2022-05-25 [2] CRAN (R 4.3.0)
# utf8                     1.2.3      2023-01-31 [2] CRAN (R 4.3.0)
# vctrs                    0.6.2      2023-04-19 [2] CRAN (R 4.3.0)
# viridis                  0.6.2      2021-10-13 [2] CRAN (R 4.3.0)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    2.5.0      2022-03-03 [2] CRAN (R 4.3.0)
# WriteXLS                 6.4.0      2022-02-24 [2] CRAN (R 4.3.0)
# xfun                     0.39       2023-04-20 [2] CRAN (R 4.3.0)
# XML                      3.99-0.14  2023-03-19 [2] CRAN (R 4.3.0)
# xml2                     1.3.4      2023-04-27 [2] CRAN (R 4.3.0)
# xtable                   1.8-4      2019-04-21 [2] CRAN (R 4.3.0)
# XVector                  0.39.0     2022-12-20 [2] Bioconductor
# yaml                     2.3.7      2023-01-23 [2] CRAN (R 4.3.0)
# yulab.utils              0.0.6      2022-12-20 [2] CRAN (R 4.3.0)
# zlibbioc                 1.45.0     2022-12-20 [2] Bioconductor
# 
# [1] /Users/lailarad/Library/R/arm64/4.3/library
# [2] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library

