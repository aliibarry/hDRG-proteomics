library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(viridis)
library(ComplexHeatmap)
library(MSFA)
library(psych)
library(multcomp)
source("./src/functions.R") #includes start_msfa_amb with more iterations

#dir.create("./output/MSFA")
PATH_results = "./output/MSFA/"

set.seed(52)

#-------------------------------------------------------------------------------

# can skip down to plotting by loading MSFA_data.RData

proteomics  <- read.csv("./data/processed/matrix-for-limma.csv", row.names = 1)
colData     <- read.csv("./data/processed/colData-for-limma.csv")

proteomics <- proteomics %>%
  distinct(genes, .keep_all = TRUE)

row.names(proteomics) <- proteomics$genes
proteomics$genes <- NULL

rna         <- read.csv("./data/multi-omics/ray2023-qnTPMs.csv", header=TRUE, check.names=FALSE, row.names = 1) #../hdrg_immune/processing/ray2023-qnTPMs.csv
rna_colData <- read.csv("./data/multi-omics/ray2023-colData.csv", header=TRUE, check.names=FALSE)

DARs <- read.csv("./data/multi-omics/DARs_Franco-Enzastiga2024.csv", header = TRUE)
DARs <- DARs$gene

#-------------------------------------------------------------------------------

# convert ensembl ids with biomaRt
rna$ensembl_gene_id <- rownames(rna)

ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_data <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters = 'ensembl_gene_id', values = rownames(rna), mart = ensembl)
gene_data <- gene_data[!duplicated(gene_data[,c('ensembl_gene_id')]),]

rna <- merge(rna, gene_data, by = "ensembl_gene_id")
rna <- rna %>% distinct(hgnc_symbol, .keep_all = TRUE)
rownames(rna) <- rna$hgnc_symbol
rna$ensembl_gene_id <- NULL
rna$hgnc_symbol <- NULL
rna$description <- NULL

rm(ensembl, gene_data)

rna        <- rna[complete.cases(rna), ]
proteomics <- proteomics[complete.cases(proteomics), ]

# Find the indices of rows with non-zero variance in both data frames
rna_var            <- apply(rna, 1, var)
proteomics_var     <- apply(proteomics, 1, var)

rna_var_idx        <- which(rna_var != 0)
proteomics_var_idx <- which(proteomics_var != 0)

rna        <- rna[rna_var_idx, ]
proteomics <- proteomics[proteomics_var_idx, ]

# MSFA requires matching columns,`P` (here, genes...will transpose matrix below).
common_rownames   <- intersect(rownames(rna), rownames(proteomics))
common_rownames   <- common_rownames[common_rownames %in% DARs] #reduce matrix size by DARs

subset_rna        <- rna[common_rownames, ]
subset_proteomics <- proteomics[common_rownames, ]

rna_matrix        <- as.matrix(subset_rna)
proteomics_matrix <- as.matrix(subset_proteomics)


#-------------------------------------------------------------------------------

# Prep for MSFA
data1 <- t(rna_matrix)
data2 <- t(proteomics_matrix)

common_genes <- intersect(colnames(data1), colnames(data2))

data1 <- data1[, common_genes]
data2 <- data2[, common_genes]

#check suitability for the factor analysis 
KMO(scale(data1)) #0.5
KMO(scale(data2)) #0.5

pca_1 <- prcomp(data1)
pca_2 <- prcomp(data2)

screeplot(pca_1)
screeplot(pca_2)

blocks <- list(
  RNA = data1,
  Proteomics = data2
)

str(blocks)
any(is.na(blocks))

#-------------------------------------------------------------------------------

# long runtime, can also just load data in next section
start_value <- start_msfa_amb(X_s = blocks, k = 3, j_s = c(4,3), robust = FALSE, constraint = "block_lower2") # initiate with more iterations, 10000 cap like De Winter et al.

mle         <-  ecm_msfa(X_s = blocks, 
                 #nIt = 1000, #debug with smaller nIt, default 50000
                 start = start_value, 
                 robust = FALSE) 

save(blocks, data1, data2, start_value, rna_colData, colData, mle, file = "./data/multi-omics/MSFA_full.RData")

#-------------------------------------------------------------------------------

load("./data/multi-omics/MSFA_full.RData")
str(mle)

# #plot estimated matrix
# gplots::heatmap.2(mle$Phi,
#                   dendrogram='row', 
#                   Rowv=TRUE, 
#                   Colv=FALSE,
#                   trace='none', 
#                   density.info="none", 
#                   col=heat.colors(256))

# Extract factor loadings and specific variances for the two datasets
lambda_s1 <- mle$Lambda_s[[1]]
lambda_s2 <- -mle$Lambda_s[[2]] #invert sign to match direction
psi_s1 <- mle$psi_s[[1]]
psi_s2 <- mle$psi_s[[2]]

# Convert specific variances to a diagonal matrix
Psi1 <- diag(psi_s1)
Psi2 <- diag(psi_s2)

# Calculate the factor scores using the bartlett method, https://online.stat.psu.edu/stat505/lesson/12/12.12
# Factor scores = (Lambda' * Psi^(-1) * Lambda)^(-1) * Lambda' * Psi^(-1) * data

# For dataset 1
lambda_s1_inv  <- solve(t(lambda_s1) %*% solve(Psi1) %*% lambda_s1)
factor_scores1 <- lambda_s1_inv %*% t(lambda_s1) %*% solve(Psi1) %*% t(scale(data1))

# For dataset 2
lambda_s2_inv  <- solve(t(lambda_s2) %*% solve(Psi2) %*% lambda_s2)
factor_scores2 <- lambda_s2_inv %*% t(lambda_s2) %*% solve(Psi2) %*% t(scale(data2))

# Transpose the factor scores to match samples x factors
factor_scores1 <- t(factor_scores1)
factor_scores2 <- t(factor_scores2)

# extract sex info from metadata
data1_meta <- data.frame(sample = rna_colData$sample_id, sex = rna_colData$sex)
data1_meta <- unique(data1_meta)
data1_meta <- data1_meta[match(rownames(data1), data1_meta$sample), ]
sex1 <- data1_meta$sex

data2_meta <- data.frame(sample = colData$mergedID, sex = colData$Sex)
data2_meta <- unique(data2_meta)
data2_meta <- data2_meta[match(rownames(data2), data2_meta$sample), ]
sex2 <- data2_meta$sex

#-------------------------------------------------------------------------------

# Create data frames for ggplot
df1 <- data.frame(Factor1 = factor_scores1[, 1], Factor2 = factor_scores1[, 2], Factor3 = factor_scores1[, 3], Sex = sex1)
df2 <- data.frame(Factor1 = factor_scores2[, 1], Factor2 = factor_scores2[, 2], Factor3 = factor_scores2[, 3], Sex = sex2)

# Plot for Dataset 1
p1 <- ggplot(df1, aes(x = Factor2, y = Factor3, color = Sex)) 
p1 <- p1 + geom_point() 
p1 <- p1 + theme_bw() + ggtitle('RNA') +
  xlab('Factor 2') +
  ylab('Factor 3')

# Plot for Dataset 2
p2 <- ggplot(df2, aes(x = Factor2, y = Factor3, color = Sex)) 
p2 <- p2 + geom_point() 
p2 <- p2 + theme_bw() + ggtitle('Proteomics') +
  xlab('Factor 2') +
  ylab('Factor 3')

pdf(paste0(PATH_results, "MSFA_scatter.pdf"), height = 3, width = 7)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

#-------------------------------------------------------------------------------

df1$Dataset <- "rna"
df2$Dataset <- "proteomics"

# Combine dataframes
df <- bind_rows(df1, df2)
df <- pivot_longer(df, cols = starts_with("Factor"), names_to = "Factor", values_to = "Value")

# Plot using ggplot2 with facet_grid
g <- ggplot(df, aes(x = Sex, y = Value)) 
g <- g + geom_boxplot(aes(fill = Sex))
g <- g + theme_bw() + facet_grid(Dataset ~ Factor) #, scales = "free"

pdf(paste0(PATH_results, "MSFA_boxplot_byfactor.pdf"), height = 4, width = 6)
print(g)
dev.off()

# df$Sex     <- as.factor(df$Sex)
# df$Dataset <- as.factor(df$Dataset)
# df$Factor  <- as.factor(df$Factor)

df$Interaction <- interaction(df$Sex, df$Factor)

# Perform ANOVA
anova_model <- aov(Value ~ Interaction, data = df) #
summary(anova_model)

contrasts <- glht(anova_model, linfct = mcp(Interaction = c(
  "M.Factor1 - F.Factor1 = 0",
  "M.Factor2 - F.Factor2 = 0",
  "M.Factor3 - F.Factor3 = 0"
)))

summary(contrasts)

anova_summary <- summary(anova_model)
anova_table   <- anova_summary[[1]]

contrast_summary <- summary(contrasts)

contrast_df <- data.frame(
  comparison = names(contrast_summary$test$coefficients),
  estimate = contrast_summary$test$coefficients,
  `std.err` = contrast_summary$test$sigma,
  `t.value` = contrast_summary$test$tstat,
  `padj` = contrast_summary$test$pvalues
)

# Write the data frame to a CSV file
write.csv(anova_table, paste0(PATH_results, "MSFA_anova.csv"))
write.csv(contrast_df, paste0(PATH_results, "MSFA_posthoc.csv"))

#-------------------------------------------------------------------------------

# GSEA on loadings for relevant factors 2 and 3
loadings <- mle$Phi
rownames(loadings) <- colnames(data1)

head(loadings)

gene_loadings <- loadings[,3] #select factor 2-3
gene_loadings <- gene_loadings[order(gene_loadings, decreasing=TRUE)]
plot(gene_loadings)

# Convert the named vector to a dataframe
gene_loadings <- data.frame(
  Gene = names(gene_loadings),
  Value = as.numeric(gene_loadings)
)

# for top/bottom ranked genes
gene_loadings <- gene_loadings %>%
  arrange(Value) %>%
  mutate(Rank = row_number())

# Identify the top and bottom 5%
n <- nrow(gene_loadings)
top_5_percent <- ceiling(n * 0.10)
bottom_5_percent <- top_5_percent

gene_loadings <- gene_loadings %>%
  mutate(Label = case_when(
    Rank <= bottom_5_percent ~ Gene,
    Rank > (n - top_5_percent) ~ Gene,
    TRUE ~ NA_character_
  ))


# Add a label column for significant pathway genes
tnfa <- c("CCNL1", "CEBPD", "PLPP3", "ABCA1", "SERPINB8")

gene_loadings <- gene_loadings %>%
  mutate(Label = case_when(
    Gene %in% tnfa ~ Gene,
    TRUE ~ NA_character_
  ))

g <- ggplot(gene_loadings, aes(x = Rank, y = Value)) 
g <- g + geom_point(color = "black") 
g <- g + ggrepel::geom_text_repel(aes(label = Label), 
                                  vjust = 2, hjust = 1, max.overlaps = 15, size = 4,
                                  box.padding = 0.35, point.padding = 0.8) 
g <- g + labs(title = "Factor 3", x = "Rank", y = "Factor Loading") 
g <- g + theme_minimal()

print(g)

#-------------------------------------------------------------------------------

library(ggrepel)
library(msigdbr)
library("clusterProfiler")
library(ggraph)

g <- ggplot(gene_loadings, aes(x = Rank, y = Value)) 
g <- g + geom_point(color = "black") 
g <- g + geom_label_repel(data = gene_loadings, mapping = aes(x = Rank, y = Value, label = Label), 
                          stat = "identity", vjust = 0.75, color = "#2c89a0ff", segment.color = 'grey50', 
                          force = 50, box.padding = 0.35, point.padding = 0.5, abel.size = 0.15, max.overlaps = 25)
g <- g + labs(title = "Factor 3", x = "Rank", y = "Factor Loading") 
g <- g + theme_minimal()
print(g)

pdf(paste0(PATH_results, "/MSFA_sharedloadings_factor3_TNFA.pdf"), height = 4, width = 5)
print(g)
dev.off()

write.csv(gene_loadings, paste0(PATH_results, "MSFA_factor3_loadings.csv"))

#-------------------------------------------------------------------------------

# GSEA on ranked loadings
hallmark_sets = msigdbr(species = "Homo sapiens", category = "H")

input <- gene_loadings

# rank
filtered_dge <- input %>%
  dplyr::arrange(dplyr::desc(abs(Value))) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

lfc_vector <- filtered_dge$Value
names(lfc_vector) <- filtered_dge$Gene

# Sort loading values in descending order
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

set.seed(52)

gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 5, 
  #maxGSSize = 500,
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(hallmark_sets, gs_name, gene_symbol)
)

gsea_result_df <- data.frame(gsea_results@result)
write.csv(gsea_result_df, file = paste0(PATH_results, "GSEA_results_factor3.csv"))

g <- ggplot(gsea_result_df, aes(x=(Description), y=NES, colour=p.adjust, size = 5)) #, size=setSize
g <- g + geom_point() + theme_bw() + ggtitle("MSFA, factor 3") #update factor
g <- g + ylim(-2, 2)
g <- g + theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
               axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
               legend.text = element_text(size=10), 
               axis.title.x = element_blank(),
               plot.title=element_text(size=rel(1), hjust = 1)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
  labs(y="Enrichment Score", colour="p value", size = "") #, size="Count"

pdf(paste0(PATH_results, "GSEA_dots_hallmark_factor3_sig.pdf"), height = 5, width = 3)
print(g)
dev.off()

#-------------------------------------------------------------------------------

key_genes <- data.frame(genes = gsea_result_df$core_enrichment, pathway = gsea_result_df$Description)
key_genes <- key_genes$genes[key_genes$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"]
key_genes <- strsplit(key_genes, "/")[1]

mycats <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

pother <- cnetplot(gsea_results, showCategory = mycats,  color_category="black", foldChange = gsea_results@geneList) 

ggraph(pother$data) +
  geom_edge_link(alpha=.4) + 
  geom_node_point(data = pother$data, aes(size= abs(color), fill = color), shape=21,  alpha = 0.9) +
  scale_size_continuous(range = c(4, 14)) +
  scale_fill_gradientn(colours = c("#1F8F89", "white", "#EE5A45"), values = scales::rescale(c(-3.5, 0, 3.5)), limits=c(-3.5, 3.5)) +
  geom_node_text(data = 
                   filter(pother$data, name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INFLAMMATORY_RESPONSE")), 
                 aes(label = name), size = 5, repel = T, max.overlaps = Inf) +
  geom_node_text(data = filter(pother$data, !name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_INFLAMMATORY_RESPONSE")),
                 aes(label = name), size = 3, repel = T, max.overlaps = Inf) +
  theme_graph(base_family = 'Helvetica') + labs(title = "", size="abs(NES)", fill="NES")


pdf(paste0(PATH_results, "GSEA_hallmark_network_factor3.pdf"), height = 5, width = 6)
ggraph(pother$data) +
  geom_edge_link(alpha=.4) + 
  geom_node_point(data = pother$data, aes(size= abs(color), fill = color), shape=21,  alpha = 0.9) +
  scale_size_continuous(range = c(4, 14)) +
  scale_fill_gradientn(colours = c("#1F8F89", "white", "#EE5A45"), values = scales::rescale(c(-1, 0, 1)), limits=c(-1, 1)) +
  geom_node_text(data = 
                   filter(pother$data, name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")), 
                 aes(label = name), size = 5, repel = T, max.overlaps = Inf) +
  geom_node_text(data = filter(pother$data, !name %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE")),
                 aes(label = name), size = 3, repel = T, max.overlaps = Inf) +
  theme_graph(base_family = 'Helvetica') + labs(title = "", size="abs(loading)", fill="Loading")
dev.off()

