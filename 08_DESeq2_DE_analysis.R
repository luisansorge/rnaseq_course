# Loading libaries used
library(DESeq2)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(writexl)

# With the help of the DESeq2 Tutorial from Ashley Valentina Schwartz and DESEQ2 R Tutorial by Lauren Ashlock.
# https://ashleyschwartz.com/posts/2023/05/deseq2-tutorial
# https://lashlock.github.io/compbio/R_presentation.html


# -----------------------------
# 5. Exploratory data analysis
# -----------------------------


# Loads modified (removed the first line and the columns containing Chr, Start, End, Strand and Length) featureCounts output.
counts <- read.table("modified_counts.txt", header = FALSE, row.names = 1)

# Creating a metadata dataframe (colData) describing the experimental conditions for each sample.
# The "sample" column uses the column names from the counts table.
# The "condition" column specifies whether the sample is from an "infected" or "control" group.
# The "tissue" column specifies the tissue type "lung" or "blood" from which each sample was taken.
colData <- data.frame(
  sample = colnames(counts),
  condition = c(rep("infected", 5), rep("control", 3), rep("infected", 5), rep("control", 3)),
  tissue = c(rep("lung", 8), rep("blood", 8))
)

# Set the row names of the metadata dataframe (colData) to match the column names of the counts matrix.
# This ensures that each sample's metadata is correctly linked to the corresponding columns in the count data.
row.names(colData) <- colnames(counts)

# Creates a DESeqDataSet object by combining the raw count data (counts) and sample metadata (colData).
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = colData, 
                              design = ~ tissue + condition)

# Runs the DESeq function on the DESeqDataSet object (dds), estimating size factors, dispersions, and fits the negative binomial model to the count data.
dds <- DESeq(dds)

# vst() removes the dependency of the variance on the mean.
# For visualisation the transformation ignores the experimental design therefore blind=TRUE.
dds_vst <- vst(dds, blind = TRUE)

# Principle component analysis (PCA) plot with default settings.
# intgroup is a vector specifying metadata columns (from colData) to colour and label the points.
plotPCA(dds_vst, intgroup = c("condition", "tissue"))


# ------------------------------------
# 6. Differential expression analysis
# ------------------------------------


# Extracts the results from the DESeqDataSet object (dds).
# results() retrieves a table of statistics, which includes log2 fold changes, p-values and adjusted p-values (padj).
res <- results(dds)
summary(res)

# Pairwise contrast of infected and control lung tissue.
# Selection of only lung sample data and updating the design formula for dds_lung to only focus on the condition comparison.
dds_lung <- dds[, colData(dds)$tissue == "lung"]
design(dds_lung) <- ~ condition
dds_lung <- DESeq(dds_lung)
res_lung <- results(dds_lung, contrast = c("condition", "infected", "control"))
summary(res_lung)

# Counts the number of genes specific to lung tissue with an adjusted p-value (padj) less than 0.05, indicating statistically significant differential expression.
sum(res_lung$padj < 0.05, na.rm=TRUE) #10287

# Count the number of genes specific to lung tissue which are up- and down-regulated.
upregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange > 0, na.rm = TRUE)
upregulated #3860
downregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange < 0, na.rm = TRUE)
downregulated #6427

# Volcano plot to visualise the results of differential expression for lung samples.
# Significant genes with adjusted p-value (padj) less than 0.05 in dodgerblue2 color.
# Genes that are both significantly differentially expressed (padj < 0.05) and have a log2FC greater than 2 or less than -2 in firebrick3 color.
par(mfrow=c(1,1))
with(res_lung, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot of differentially expressed genes\nin Infected vs Control lung tissue"))
with(subset(res_lung, padj < 0.05), points(log2FoldChange, -log10(pvalue), pch=20, col="dodgerblue2"))
with(subset(res_lung, padj < 0.05 & abs(log2FoldChange) > 2), points(log2FoldChange, -log10(pvalue), pch=20, col="firebrick3"))

# Adding of an additional column that identifies a gene as up-regulated, down-regulated or unchanged.
data <- data.frame(res_lung)
data <- data %>%
  mutate(
    Expression = case_when(log2FoldChange >= log(1) & padj <= 0.05 ~ "Up-regulated",
                           log2FoldChange <= -log(1) & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

# Filtering the column Up-regulated and Down-regulated and sorting by the adjusted p-value (padj) and getting the top 10 most up- and down-regulated genes. 
top <- 10
top_genes <- bind_rows(
  data %>%
    filter(Expression == 'Up-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top),
  data %>%
    filter(Expression == 'Down-regulated') %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    head(top)
)

# Creates a data frame just holding the top 10 up-regulated genes.
Top_Hits = head(arrange(data,pvalue),10)

# Create a volcano plot using ggplot2 to visualise the differential expression results in lung tissue.
# Specifically colouring up-regulated genes in dodgerblue2, black for genes that are unchanged, and firebrick3 for down-regulated genes.
ggvolcano <- ggplot(data, aes(log2FoldChange, -log(pvalue,10))) + # -log10 conversion
  geom_point(aes(color = Expression), size = 3/5) +
  xlab(expression("log"[2]*"FoldChange")) +
  ylab(expression("-log"[10]*"P-Value")) +
  scale_color_manual(values = c("dodgerblue2", "black", "firebrick3")) +
  xlim(-20, 10) +
  ylim(0, 300) 
ggvolcano

# Creation of new factor variable "group" that combines the tissue type and condition.
dds$group <- factor(paste0(dds$tissue, dds$condition))

# Updates the design formula to model the differential expression based on the "group" defined above.
design(dds) <- ~ group

# Runs the DESeq function on the updated DESeqDataSet object (dds), estimating size factors, dispersions, and fits the negative binomial model to the count data.
dds <- DESeq(dds)

# Orders to show lowest adjusted p-values first (most statistically significant differential expression).
# Ranking all genes based on their significance, regardless of the tissue type or condition.
res_pval <- res[order(res$padj),]
head(res_pval)

# Ranking genes for just the lung tissue based on their significance.
res_lung_pval <- res_lung[order(res_lung$padj),]
head(res_lung_pval)

# Relevant genes with the lowest adjusted p-value (padj).
par(mfrow=c(2,2))
plotCounts(dds, gene="ENSMUSG00000028270", intgroup="group", main = expression(italic("Gbp2")))
plotCounts(dds, gene="ENSMUSG00000105504", intgroup="group", main = expression(italic("Gbp5")))
plotCounts(dds, gene="ENSMUSG00000052776", intgroup="group", main = expression(italic("Oas1a")))
plotCounts(dds, gene="ENSMUSG00000078922", intgroup="group", main = expression(italic("Tgtp1")))

# Relevant genes based on the paper by Singhania et al. (2019).
par(mfrow=c(2,3))
#plotCounts(dds, gene="ENSMUSG00000048806", intgroup="group", main = "Ifnb1")
#plotCounts(dds, gene="ENSMUSG00000034459", intgroup="group", main = "Ifit1")
plotCounts(dds, gene="ENSMUSG00000074896", intgroup="group", main = "Ifit3")
plotCounts(dds, gene="ENSMUSG00000015947", intgroup="group", main = "Fcgr1")
#plotCounts(dds, gene="ENSMUSG00000032690", intgroup="group", main = "Oas2")
#plotCounts(dds, gene="ENSMUSG00000032661", intgroup="group", main = "Oas3")
plotCounts(dds, gene="ENSMUSG00000041827", intgroup="group", main = "Oasl1")
#plotCounts(dds, gene="ENSMUSG00000029561", intgroup="group", main = "Oasl2")
plotCounts(dds, gene="ENSMUSG00000000386", intgroup="group", main = "Mx1")
plotCounts(dds, gene="ENSMUSG00000040033", intgroup="group", main = "Stat2")
#plotCounts(dds, gene="ENSMUSG00000025498", intgroup="group", main = "Irf7")
plotCounts(dds, gene="ENSMUSG00000002325", intgroup="group", main = "Irf9")

# Additional genes.
#plotCounts(dds, gene="ENSMUSG00000029417", intgroup="group", main = "Cxcl9")
#plotCounts(dds, gene="ENSMUSG00000068606", intgroup="group", main = "Gm4841")
#plotCounts(dds, gene="ENSMUSG00000039699", intgroup="group", main = "Batf2")
#plotCounts(dds, gene="ENSMUSG00000016496", intgroup="group", main = "Cd274")


# -------------------------------
# 7. Overrepresentation analysis
# -------------------------------


# Perform differential expression analysis between the infected and control groups for lung tissue.
# "contrast" specifies which conditions are being compared.
res_group_lung <- results(dds, contrast = c("group", "lunginfected", "lungcontrol"))
summary(res_group_lung)

# Perform differential expression analysis between the infected and control groups for blood.
res_group_blood <- results(dds, contrast = c("group", "bloodinfected", "bloodcontrol"))
summary(res_group_blood)

# Extraction of gene ids of genes that are significantly differentially expressed for the overall dataset (res).
geneids <- rownames(res[res$padj < 0.05 & !is.na(res$padj),])

# Extraction of gene ids of genes that are significantly differentially expressed specifically for lung tissue comparisons.
geneidslung <- rownames(res_group_lung[res_group_lung$padj < 0.05 & !is.na(res_group_lung$padj),])

# Extraction of gene ids of genes that are significantly differentially expressed specifically for blood comparisons.
geneidsblood <- rownames(res_group_blood[res_group_blood$padj < 0.05 & !is.na(res_group_blood$padj),])

# Contains all gene ids from the overall results, regardless of significance.
all.geneids <- rownames(res)

# Perform Gene Ontology (GO) overrepresentation analysis on the differentially expressed genes from the lung tissue (geneidslung).
ego_lung <- enrichGO(
  gene = geneidslung,            # gene ids for genes differentially experessed in lung tissue
  universe = all.geneids,        # all gene ids (background)
  OrgDb = org.Mm.eg.db,          # mouse annotation database
  ont = "BP",                    # subontology biological process
  keyType = "ENSEMBL",           # gene id format
  pAdjustMethod = "BH",          # Adjust p-values using Benjamini-Hochberg method
  pvalueCutoff = 0.05,           # Significant GO terms
  qvalueCutoff = 0.2,            # adjusted q-value cutoff
  readable = TRUE                # converts ENSEMBL ids to gene symbols
)

# Perform Gene Ontology (GO) overrepresentation analysis on the differentially expressed genes from the blood (geneidsblood).
ego_blood <- enrichGO(
  gene = geneidsblood,           # gene ids for genes differentially experessed in blood
  universe = all.geneids,        # all gene ids (background)
  OrgDb = org.Mm.eg.db,          # mouse annotation database
  ont = "BP",                    # subontology biological process
  keyType = "ENSEMBL",           # gene id format
  pAdjustMethod = "BH",          # Adjust p-values using Benjamini-Hochberg method
  pvalueCutoff = 0.05,           # Significant GO terms
  qvalueCutoff = 0.2,            # adjusted q-value cutoff
  readable = TRUE                # converts ENSEMBL ids to gene symbols
)

# Creates barplot and dotplot of the Top 10 GO terms in lung tissue.
barplot(ego_lung, showCategory = 10)  
dotplot(ego_lung, showCategory = 10)

# Creates barplot and dotplot of the Top 10 GO terms in blood.
barplot(ego_blood, showCategory = 10)
dotplot(ego_blood, showCategory = 10)

go_results_ego_lung <- ego_lung@result  #Lung tissue results
go_results_ego_blood <- ego_blood@result  #Blood tissue results

#Creating the results in a dataframe
go_results_ego_lung_df <- as.data.frame(go_results_ego_lung)
go_results_ego_blood_df <- as.data.frame(go_results_ego_blood)

#Creating list of results
results_list <- list(
  "GO_Lung_tissue" = go_results_ego_lung_df,
  "GO_Blood_tissue" = go_results_ego_blood_df
)

#Write all results to an Excel file
write_xlsx(results_list, path = "GO_Enrichment_Results_Lung_Blood.xlsx")
