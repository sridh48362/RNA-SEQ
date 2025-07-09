library(readr)
library(dplyr)
df <- read_tsv("combine.tsv")

df_clean <- df %>% filter(grepl("^ENSG", Features))

counts <- df_clean %>% column_to_rownames("Features")

sample_names <- colnames(counts)

conditions <- c(rep("BT", 3), 
                rep("STM", 2),  
                rep("STM_BT", 3),   
                "STM_C",            
                rep("UT", 3)) 

coldata <- data.frame(row.names = sample_names,
                      condition = factor(conditions))
library(edgeR)
keep <- rowSums(cpm(counts) > 1) >= 2
counts_filtered <- counts[keep, ]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")

res_UT_BT <- results(dds, contrast = c("condition", "UT", "BT"))
res_UT_BT <- res_UT_BT[order(res_UT_BT$padj), ]
head(res_UT_BT)

library(EnhancedVolcano)

EnhancedVolcano(res_UT_BT,
                lab = rownames(res_UT_BT),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)

library(org.Hs.eg.db)
library(AnnotationDbi)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(res_UT_BT),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add gene symbols to the DESeq2 result
res_UT_BT$gene_symbol <- gene_symbols


res_df <- as.data.frame(res_UT_BT)

# filter upregulated and downregulated genes
upregulated <- res_df %>%
  filter(log2FoldChange > 1, padj < 0.05)

downregulated <- res_df %>%
  filter(log2FoldChange < -1, padj < 0.05)


library(clusterProfiler)

# GO enrichment analysis
enriched_up <- enrichGO(gene = up_genes,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)

enriched_down <- enrichGO(gene = down_genes,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

# Plot 
barplot(enriched_up, showCategory = 10)
barplot(enriched_down, showCategory = 10)





