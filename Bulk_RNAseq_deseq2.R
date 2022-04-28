library(tximport)
library(readr)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(apeglm)
library("BiocParallel")
register(MulticoreParam(4))
library("pheatmap")
library("tidyverse")
library(edgeR)
library(RColorBrewer)

#Set directory to location of TSV values
setwd("~/Kallisto_quant_files")
files = list.files()

#Reorder files
files = files[c(12:14,1:11)]

#Get sample names
samples=gsub("_abundance.tsv","",files)
names(files) = samples

#Make key to convert ensembl transcript ids to gene names
edb = EnsDb.Hsapiens.v86
tx = transcripts(edb, columns = c("tx_id", "gene_name"), return.type = "DataFrame")

#Import transcript abundance data
txi = tximport(files, type = "kallisto", tx2gene = tx, ignoreTxVersion = TRUE)

#Set up metadata
cell_type = as.factor(c(rep("artery",2),
                        rep("vein", 2)))

#Add metadata columns 
comparisons = c("artery", "vein")
meta_data = data.frame(sample = samples, cell_type = cell_type)
meta_data[comparisons] = 0 
for (i in 1:length(comparisons)){
  meta_data[which(meta_data$cell_type==comparisons[i]),comparisons[i]]=1
}

meta_data[comparisons]=lapply(meta_data[comparisons], factor)
class(meta_data)

#General dds
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta_data,
                                design = ~ cell_type)
as.data.frame(colData(dds))
dds<- DESeq(dds)
class(dds)
rld<-rlog(dds)

#Make DESeqDataSet for artery comparison
dds_artery = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ artery)

#Run differential expression analysis with default parameters
dds_artery = DESeq(dds_artery)
resultsNames(dds_artery)
res_artery = results(dds_artery, contrast = c("artery", "1", "0"), name="artery_1_vs_0")
res_artery
res_artery_Sig <- subset(res_artery, padj < 0.01)

res_artery_df <- as.data.frame(res_artery_Sig) 
res_artery_tibble <- as.tibble(res_artery_df, rownames = "geneID") 
res_artery_tibble<-subset(res_artery_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_artery_tibble), file="DESEQ2_top_artery_p<0.01_tibble.csv" )
res_artery_top <- res_artery_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_artery_top
head(res_artery_top, 15)
res_artery_top80 <- head(res_artery_top, 80)

#Make DESeqDataSet for vein comparison 
dds_vein = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ vein)

#Run differential expression analysis with default parameters
dds_vein = DESeq(dds_vein)
res_vein = results(dds_vein, contrast = c("vein", "1", "0"))
res_vein
res_vein_Sig <- subset(res_vein, padj < 0.01)

res_vein_df <- as.data.frame(res_vein_Sig) 
res_vein_tibble <- as.tibble(res_vein_df, rownames = "geneID") 
res_vein_tibble<-subset(res_vein_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_vein_tibble), file="DESEQ2_top_vein_p<0.01_tibble.csv" )
res_vein_top <- res_vein_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_vein_top
head(res_vein_top, 15)
res_vein_top80 <- head(res_vein_top,80)

#Plot heatmap for artery vein comparison---- 

rownames<- c(res_artery_top$geneID, res_vein_top$geneID)
dds
colData(dds)
ddsColl <- collapseReplicates(dds, dds$cell_type, renameCols=FALSE)
ddsColl

colData(ddsColl)
ddsColl <- ddsColl[,order(ddsColl$cell_type)]
rldColl<-rlog(ddsColl)

pheatmap(assay(rld)[rownames,], scale="row", fontsize = 20, cluster_cols = FALSE, 
         cluster_rows = FALSE, col= colorRampPalette( rev(brewer.pal(12, "RdBu")) )(255))
pheatmap(assay(rldColl)[rownames,],
         scale="row", fontsize = 20, cluster_cols = F, show_rownames = F,
         border_color = "grey60", cluster_rows = FALSE, col= colorRampPalette(colors = c("#3b3bff", "#000014", "#ff0000"))(255))
pheatmap(assay(rld)[rownames,], 
         scale="row", fontsize = 20, cluster_cols = F, show_rownames = F,
         border_color = "grey60", cluster_rows = FALSE, col= colorRampPalette(colors = c("#3b3bff", "#000014", "#ff0000"))(255))
pheatmap(assay(rld)[rownames,], scale="row", fontsize = 20, cluster_cols = FALSE, show_rownames = F,
         cluster_rows = FALSE, col= colorRampPalette(colors = c("#1E90FF", "#000014", "#ff0000"))(255))

ggsave("heatmap_top50to80.tiff", plot=last_plot(), height = 7, width =5)



