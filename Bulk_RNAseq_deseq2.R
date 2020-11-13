library(tximport)
library(readr)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
BiocManager::install("apeglm")
library(apeglm)
library("BiocParallel")
register(MulticoreParam(4))
library("pheatmap")
library("tidyverse")
library(edgeR)
library(RColorBrewer)

#set directory to location of TSV values----
setwd("/Users/LTANG/Downloads/Kallisto_quant_files")
files = list.files()

#reorder files
files = files[c(12:14,1:11)]

#get sample names
samples=gsub("_abundance.tsv","",files)
names(files) = samples

#make key to convert ensembl transcript ids to gene names
edb = EnsDb.Hsapiens.v86
tx = transcripts(edb, columns = c("tx_id", "gene_name"), return.type = "DataFrame")

#import transcript abundance data
txi = tximport(files, type = "kallisto", tx2gene = tx, ignoreTxVersion = TRUE)

#set up metadata----
cell_type = as.factor(c(rep("hesc",3), rep("ps",3), rep("dlm",2), rep("artery",2),
                        rep("pre_vein",2), rep("vein",2)))

#add metadata columns where a given cell type has a 1 and all other cell types have a 0
comparisons = c("hesc", "ps", "dlm", "artery", "pre_vein", "vein")
meta_data = data.frame(sample = samples, cell_type = cell_type)
meta_data[comparisons] = 0 #[] is similar to $ to pick columns
for (i in 1:length(comparisons)){
  meta_data[which(meta_data$cell_type==comparisons[i]),comparisons[i]]=1
}

meta_data[comparisons]=lapply(meta_data[comparisons], factor)
class(meta_data)
#general dds----
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta_data,
                                design = ~ cell_type)
as.data.frame(colData(dds))
dds<- DESeq(dds)
class(dds)
rld<-rlog(dds)

#make DESeqDataSet for artery comparison----
dds_artery = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ artery)

#run differential expression analysis with default parameters
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

#make DESeqDataSet for vein comparison----
dds_vein = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ vein)

#run differential expression analysis with default parameters
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

#make DESeqDataSet for pre_vein comparison----
dds_pre_vein = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ pre_vein)

#run differential expression analysis with default parameters
dds_pre_vein = DESeq(dds_pre_vein)
res_pre_vein = results(dds_pre_vein, contrast = c("pre_vein", "1", "0"))
res_pre_vein
res_pre_vein_Sig <- subset(res_pre_vein, padj < 0.01)
res_pre_vein_df <- as.data.frame(res_pre_vein_Sig) 
res_pre_vein_tibble <- as.tibble(res_pre_vein_df, rownames = "geneID") 
res_pre_vein_tibble<-subset(res_pre_vein_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_pre_vein_tibble), file="DESEQ2_top_pre_vein_p<0.01_tibble.csv" )
res_pre_vein_top <- res_pre_vein_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_pre_vein_top
head(res_pre_vein_top, 15)
res_pre_vein_top50 <- head(res_pre_vein_top,50)

#make DESeqDataSet for dlm comparison----
dds_dlm = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ dlm)

#run differential expression analysis with default parameters
dds_dlm = DESeq(dds_dlm)
res_dlm = results(dds_dlm, contrast = c("dlm", "1", "0"))
res_dlm
res_dlm_Sig <- subset(res_dlm, padj < 0.01)
res_dlm_df <- as.data.frame(res_dlm_Sig) 
res_dlm_tibble <- as.tibble(res_dlm_df, rownames = "geneID") 
res_dlm_tibble<-subset(res_dlm_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_dlm_tibble), file="DESEQ2_top_dlm_p<0.01_tibble.csv" )
res_dlm_top <- res_dlm_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_dlm_top
head(res_dlm_top, 15)
res_dlm_top80 <- head(res_dlm_top,80)

#make DESeqDataSet for ps comparison----
dds_ps = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ ps)

#run differential expression analysis with default parameters
dds_ps = DESeq(dds_ps)
res_ps = results(dds_ps, contrast = c("ps", "1", "0"))
res_ps

res_ps_Sig <- subset(res_ps, padj < 0.01)
res_ps_df <- as.data.frame(res_ps_Sig) 
res_ps_tibble <- as.tibble(res_ps_df, rownames = "geneID")
res_ps_tibble<-subset(res_ps_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_ps_tibble), file="DESEQ2_top_ps_p<0.01_tibble.csv" )
res_ps_top <- res_ps_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_ps_top
head(res_ps_top, 15)
res_ps_top80 <- head(res_ps_top,80)

#make DESeqDataSet for hesc comparison----
dds_hesc = DESeqDataSetFromTximport(txi, colData = meta_data, design = ~ hesc)

#run differential expression analysis with default parameters
dds_hesc = DESeq(dds_hesc)
res_hesc = results(dds_hesc, contrast = c("hesc", "1", "0"))
res_hesc
res_hesc_Sig <- subset(res_hesc, padj < 0.01)
res_hesc_df <- as.data.frame(res_hesc_Sig) 
res_hesc_tibble <- as.tibble(res_hesc_df, rownames = "geneID")
res_hesc_tibble<-subset(res_hesc_tibble, log2FoldChange>1.5)
write.csv( as.data.frame(res_hesc_tibble), file="DESEQ2_top_hesc_p<0.01_tibble.csv" )

res_hesc_top <- res_hesc_tibble %>% 
  dplyr::arrange(-log2FoldChange) %>% #to ensure no conflict in function
  dplyr::select(geneID)
res_hesc_top
head(res_hesc_top, 20)
class(res_hesc_top)
write.csv( as.data.frame(res_hesc_top), file="DESEQ2_top_hesc_p<0.01.csv" )
res_hesc_top80 <- head(res_hesc_top,80)

#plot heatmap for artery comparison---- 

rownames<- c(res_hesc_top$geneID, res_ps_top$geneID, res_dlm_top$geneID, res_artery_top$geneID, res_pre_vein_top$geneID, res_vein_top$geneID)
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



