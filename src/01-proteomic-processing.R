#library(readxl)
#library(pheatmap)

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(viridis)
library(ComplexHeatmap)
library(matrixStats)
library(gridExtra)
library(stringr)

library(variancePartition)
#library(limma)

################################################################################

# library(diann)
# 
# dir.create("./data/")
# 
# # Load DIA-NN output matrix
# df <- diann_load("./report.tsv") #/path/to/diann/output
# 
# precursors <- diann_matrix(df, q = 0.01)
# 
# write.csv(precursors, "./data/precursors_all_FDR0.01.csv", 
#           na="NA",eol = "\n", row.names = T)
# 
# # calculate maxLFQ
# peptides.maxlfq <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$PG.Q.Value <= 0.01,], 
#                                 sample.header = "Run",
#                                 group.header="Stripped.Sequence", 
#                                 id.header = "Precursor.Id", 
#                                 quantity.header = "Precursor.Normalised")
# 
# peptides.maxlfq <- as.data.frame(peptides.maxlfq)
# 
# peptides.maxlfq$ModifiedSequence <- rownames(peptides.maxlfq)
# rownames(peptides.maxlfq) <- NULL
# 
# peptides.maxlfq <- peptides.maxlfq %>%
#   select(ModifiedSequence, everything())
# 
# head(peptides.maxlfq)
# 
# write.csv(peptides.maxlfq, "./data/quant_peptides.csv", 
#           na="NA", eol = "\n", row.names = FALSE)
# 
# gene.groups <- diann_maxlfq(df[df$Q.Value <= 0.01 & df$GG.Q.Value <= 0.01,], 
#                             sample.header = "Run",
#                             group.header="Genes", 
#                             id.header = "Precursor.Id", 
#                             quantity.header = "Precursor.Normalised")
# 
# write.csv(gene.groups, "./data/quant_genegroups.csv", 
#           na="NA", eol = "\n", row.names = T)

################################################################################

dir.create("./output/")
PATH_results = "./output/"

precursors <- read.csv("./data/precursors_all_FDR0.01.csv", header = TRUE, 
                       sep = ",", check.names = FALSE, row.names = 1)

peptides.maxlfq <- read.csv("./data/quant_peptides.csv", na="NA", header = TRUE, 
                        sep = ",", check.names = FALSE, row.names = 1)

gene.groups <- read.csv("./data/quant_genegroups.csv", na="NA", header = TRUE, 
                        sep = ",", check.names = FALSE, row.names = 1)

colData <- read.csv("./data/metadata.csv", sep = ",", header = TRUE)

################################################################################

# basic QC

################################################################################

pre_counts  <- colSums(!is.na(precursors))
pre_counts  <- as.data.frame(pre_counts)

sorted_indices <- order(-pre_counts$pre_counts)
sorted_counts  <- pre_counts$pre_counts[sorted_indices]
sorted_names   <- rownames(pre_counts)[sorted_indices]

pdf(file = paste(PATH_results, "counts-precursor.pdf", sep=""), width = 8, height = 6)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        ylab = "Precursor Count", 
        las = 3,
        cex.names = 0.8)
dev.off()

# GG counts
gg_counts  <- colSums(!is.na(gene.groups))
gg_counts  <- as.data.frame(gg_counts)

sorted_indices <- order(-gg_counts$gg_counts)
sorted_counts <- gg_counts$gg_counts[sorted_indices]
sorted_names <- rownames(gg_counts)[sorted_indices]

pdf(file = paste(PATH_results, "counts-gg.pdf", sep=""), width = 8, height = 8)
par(mfrow=c(2,1))
barplot(sorted_counts, 
        names.arg = sorted_names,
        #xlab = "Sample", 
        ylab = "GG Count", 
        #main = "Counts by Sample",
        las = 2,
        cex.names = 0.8)
dev.off()

pdf(file = paste(PATH_results, "counts-gg-boxplot.pdf", sep=""), width = 4, height = 4)
boxplot(sorted_counts)
dev.off()

# sample correlations
mat <- gene.groups[complete.cases(gene.groups), ]
correlation_matrix <- cor(mat)

melted_corr_matrix <- melt(correlation_matrix)

pdf(file = paste(PATH_results, "correlation.pdf", sep=""), width = 8, height = 8)
ComplexHeatmap::pheatmap(correlation_matrix,
                         main = "Correlation matrix",
                         border_color = "black",
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         col=viridis(100),
                         fontsize_row = 8,
                         fontsize_col = 8
)

dev.off()

################################################################################

# Clustering

################################################################################

data <- gene.groups

# confirm samples match colData
mat <- as.matrix(data[, which(colnames(data) %in% colData$sampleID)])
mat <- mat[complete.cases(mat), ]

rv     <- matrixStats::rowVars(mat) # calculate variance per row (ie. per gene)
select <- order(rv, decreasing=TRUE)[seq_len(min(10000, length(rv)))]

#pca <- prcomp(t(mat[select,]), center = TRUE, scale. = TRUE)
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

tissue <- as.factor(colData$Tissue)
sex    <- as.factor(colData$Sex)
age    <- as.factor(colData$Age)
replicate <- as.factor(colData$Replicate)
donor  <- as.factor(colData$donorID)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(tissue),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               labels = NULL,
               point.size = 5,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)

g2 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex),
               #ellipse = TRUE,
               ellipse.prob = 0.95,
               point.size = 5,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE, circle.prob = 0.5,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g2 <- g2 + theme(legend.position = 'right', aspect.ratio=1)

pdf(file = paste(PATH_results, "pca.pdf", sep=""), width = 8, height = 4)
g1
#grid.arrange(g1, g2, ncol = 2)
dev.off()

################################################################################

form_full <- ~ Age + Replicate + (1|Tissue) + (1|Sex) + (1|donorID)

# extract matching sample names 
mat <- as.matrix(data[, which(names(data) %in% colData$sampleID)])
rv  <- matrixStats::rowVars(mat)

select <- order(rv, decreasing=TRUE)[seq_len(min(5000, length(rv)))]

filtered_mat <- mat[select,]

varPart <- fitExtractVarPartModel(filtered_mat, form_full, colData) #update form as needed

vp <- sortCols(varPart)

g <- plotPercentBars( vp[1:10,], col = magma(11)) + 
  theme(axis.text.y = element_text(size= 10), 
        axis.title.y = element_text(size=10), axis.title.x = element_text(size= 10), 
        axis.text.x = element_text(size= 10), legend.title=element_text(size=10), 
        legend.text=element_text(size=10), plot.title=element_text(size=10, hjust = 0.5)) + 
  ggtitle("Variance Partition")

print(g)

g <- plotVarPart(sortCols(varPart), label.angle=60 )

pdf(file = paste0(PATH_results, "/variance-partitioning.pdf"))
print(g)
dev.off()

################################################################################

neurons <- read.csv("./data/neuronal-enrichment-rankings.csv", header = FALSE)
neurons <- neurons$V1

neurons <- trimws(as.character(neurons))

data <- as.data.frame(gene.groups)
data$genes <- rownames(data)

# data <- data %>%
#   separate(genes, into = paste0("genes"), sep = ";", remove = TRUE) 

data$genes <- trimws(as.character(data$genes))

data <- data %>%
  mutate(genes = gsub("\\n", "", genes)) %>%
  distinct(genes, .keep_all = TRUE)

# data$genes <- NULL
# expression_data <- data[complete.cases(data),]

expression_data <- data[data$genes %in% neurons, ]
expression_data$genes <- NULL

expression_data <- as.matrix(expression_data)
rownames(expression_data) <- data$genes[match(rownames(expression_data), rownames(data))]

scaled_expression <- t(scale(t(expression_data), center = TRUE))

ComplexHeatmap::pheatmap(expression_data,
         main = "Neuronal enrichment",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         col=viridis(100),
         fontsize_row = 8,
         fontsize_col = 8
)

match_index <- match(colnames(scaled_expression), colData$sampleID)
metadata_reordered <- colData[match_index, ]

tissue_list <- as.factor(metadata_reordered$Tissue)

# Create a color mapping for metadata
# tissue_colors <- c("F" = "grey",
#                    "M" = "darkgrey")

tissue_colors <- c("nerve.root" = "grey",
                   "ganglia" = "#712b4a")

# Assign colors to metadata levels
col_fun <- tissue_colors[tissue_list]

scaled_expression <- t(scale(t(expression_data)))

# Create Heatmap
ht_list <- Heatmap(scaled_expression,
                   #name = "Expression",
                   col=viridis(100),
                   clustering_distance_columns = "manhattan",
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = FALSE, #set to TRUE to double check colour legend
                   row_title = "Genes",
                   row_dend_side = "left",
                   top_annotation = HeatmapAnnotation(tissue = tissue_list, col = list(tissue = col_fun))
                   )

draw(ht_list, heatmap_legend_side = "right")

pdf(file = paste0(PATH_results, "/neuronal-heatmap.pdf"), height = 5, width = 7)
draw(ht_list, heatmap_legend_side = "right")
dev.off()


################################################################################

colData$mergedID <- paste0(colData$donorID, sep=".", colData$Tissue)

data_trans  <- as.data.frame(log2(gene.groups))
ranked_data <- data_trans %>% mutate_all(rank)

merged_data <- data_trans %>%
  tibble::rownames_to_column(var = "proteins") %>%
  pivot_longer(cols = -proteins, names_to = "sampleID", values_to = "Log2Intensity") %>%
  left_join(ranked_data %>% 
              tibble::rownames_to_column(var = "proteins") %>%
              pivot_longer(cols = -proteins, names_to = "sampleID", values_to = "Rank"),
            by = c("proteins", "sampleID")) %>%
  left_join(colData, by = "sampleID") %>%
  separate(proteins, into = "genes", sep = ";", remove = FALSE) %>%
  mutate(genes = trimws(as.character(genes)))

df <- merged_data %>%
  group_by(proteins, genes, mergedID) %>%
  summarize(expression = mean(Log2Intensity, na.rm = TRUE))

df <- df %>% pivot_wider(names_from = mergedID, values_from = expression)
df <- as.data.frame(df)

df <- df %>%
  mutate(proteins = gsub("\\n", "", proteins)) %>%
  distinct(proteins, .keep_all = TRUE)

head(df)

# reduce colData for merge
to.keep <- c("donorID", "Sex", "DRG", "Age", "Tissue", "mergedID")
df_meta <- colData
df_meta <- colData[, colnames(colData) %in% to.keep]
df_meta <- df_meta[!duplicated(df_meta), ]

write.csv(df,      "./data/processed/matrix-for-limma.csv", na="NA", eol = "\n", row.names = FALSE)
write.csv(df_meta, "./data/processed/colData-for-limma.csv", na="NA", eol = "\n", row.names = FALSE)

################################################################################

ranked_mean <- merged_data %>%
  select(proteins, genes, Tissue, Log2Intensity) %>%
  group_by(proteins, Tissue, genes) %>%
  summarize(mean_Log2Intensity = mean(Log2Intensity, na.rm = TRUE), .groups = 'drop') %>%
  group_by(Tissue) %>%
  mutate(Rank = rank(-mean_Log2Intensity, ties.method = "first")) %>%
  ungroup()

# superimpose key genes on plot
named_mean <- ranked_mean %>%
  filter(genes %in% neurons) %>%
  group_by(genes, Tissue) %>%
  filter(Rank == min(Rank)) %>%
  ungroup()

head(named_mean)

# dynamic range by group
g <- ggplot() +
  geom_line(data = ranked_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), color = "#FFE6BF", linewidth = 2) +
  geom_point(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity), size = 2, shape = 21, fill = "#21130d", colour = "#21130d")
g <- g + ggrepel::geom_label_repel(data = named_mean, mapping = aes(x = Rank, y = mean_Log2Intensity, label = genes), 
                                   stat = "identity", color = "#21130d", segment.color = 'grey50', 
                                   force = 50, box.padding = 0.35, point.padding = 0.5, max.overlaps = 45) 
g <- g + facet_grid(Tissue ~ .) + theme_bw() 
g <- g + labs(title = "Key proteins", x = "Rank", y = "Log2 Intensity")

print(g)

pdf(file = paste(PATH_results, "dynamic-range.pdf", sep=""), width = 6, height = 5)
print(g)
dev.off()

################################################################################






















# sPCA based on DEGs and equiv for sex in hDRG

DARs <- read.csv("./data/multi-omics/DARs_Franco-Enzastiga2024.csv")
filtered_DARs <- DARs %>% filter(!(Chr %in% c("chrX", "chrY")))

head(filtered_DARs)
head(mean_mat)

# reorder colData to match matrix
index   <- match(colnames(mean_mat), colData$mergedID)
colData <- colData[index, ]
sex     <- as.factor(colData$Sex) #sex factor as needed

mat <- mean_mat[rownames(mean_mat) %in% filtered_DARs$gene, ] #filtered_
mat <- mat[complete.cases(mat), ]

pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex),
               point.size = 5,
               ellipse = FALSE,
               ellipse.prob = 0.95,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)
g1 <- g1 + ggtitle("sPCA on ATAC-seq DARs")

print(g1)

pdf(file = paste(PATH_results, "sPCA_DARs.pdf", sep=""), width = 4, height = 4)
print(g1)
dev.off()

# extract eigengenes from this to see what drives effect
eigens <- pca$x[,1] #select PC1 frm prcomp output

eigengenes <- pca$rotation[,1]
eigengenes <- eigengenes[order(eigengenes, decreasing=TRUE)]
plot(eigengenes)

# Convert the named vector to a dataframe
eigengenes <- data.frame(
  Gene = names(eigengenes),
  Value = as.numeric(eigengenes)
)

eigengenes <- eigengenes %>%
  arrange(Value) %>%
  mutate(Rank = row_number())

# Identify the top and bottom 5%
n <- nrow(eigengenes)
top_5_percent <- ceiling(n * 0.05)
bottom_5_percent <- top_5_percent

# Add a label column
eigengenes <- eigengenes %>%
  mutate(Label = case_when(
    Rank <= bottom_5_percent ~ Gene,
    Rank > (n - top_5_percent) ~ Gene,
    TRUE ~ NA_character_
  ))

g <- ggplot(eigengenes, aes(x = Rank, y = Value)) +
  geom_point(color = "black") 
g <- g + ggrepel::geom_text_repel(aes(label = Label), 
                                  vjust = -0.5, hjust = 1, max.overlaps = 25, size = 4,
                                  box.padding = 0.35, point.padding = 0.5) 
g <- g + labs(title = "Eigengenes",
       x = "Rank",
       y = "Eigengene Value") 
g <- g + theme_minimal()

print(g)

pdf(paste0(PATH_results, "/ModuleEigengene_DARs_ranked.pdf"), height = 4, width = 4)
print(g)
dev.off()

# sex comparison
eigengenes <- pca$rotation[,1]
eigengenes <- eigengenes[order(eigengenes, decreasing=TRUE)]

input <- as.data.frame(eigengenes)
input$gene <- rownames(input)
input <- input[c("gene", "eigengenes")]

# ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
# gene_data <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description'), filters = 'ensembl_gene_id', values = rownames(input), mart = ensembl)
# gene_data <- gene_data[!duplicated(gene_data[,c('ensembl_gene_id')]),]
# rownames(gene_data) <- gene_data$ensembl_gene_id
# input$symbol <- gene_data[rownames(input),]$hgnc_symbol
# input$description <- gene_data[rownames(input),]$description

write.csv(input, file = paste0(PATH_results, "/DARs_PC1eigens.csv"))

df <- data.frame(eigen=eigens,
                 Sex = colData$Sex)

teigens <- df %>% group_by(Sex, eigen)

bartlett.test(eigen ~ Sex, data = teigens) #test for equal variance (fails, use Welch's, not aov)

df <- as.data.frame(teigens)

# rstatix::cohens_d(data = df, formula = eigen ~ Sex, var.equal = FALSE)
# 
# stat.test <- df %>%
#   rstatix::t_test(eigen ~ Sex) %>%
#   rstatix::add_significance()
# stat.test

oneway.test(eigen ~ Sex, data = df, var.equal = TRUE)

# BH <- pairwise.t.test(teigens4W$eigen, teigens4W$Pop_Cond,
#                 p.adjust.method = "BH", pool.sd = FALSE)
# 
# pairwise.wilcox.test(teigens$eigen, teigens$Sex,
#                      p.adjust.method = "BH")
# 
# kruskal.test(eigen ~ Sex, data = teigens)

##### change data for plotting (3d v 4w)
g <- ggplot(data=teigens, aes(x=interaction(Sex), y=eigen)) 
g <- g + geom_boxplot(coef = 1.5, aes(fill=Sex)) 
# g <- g + scale_fill_viridis_d() #option = "mako"
g <- g + theme_bw() + theme(
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.x = element_text(size=10, angle = 45, hjust= 1),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank())
g <- g + geom_hline(yintercept = 0)

print(g)

pdf(paste0(PATH_results, "/ModuleEigengene_DARs_boxplot.pdf"), height = 3, width = 3)
print(g)
dev.off()

################################################################################

DEGs <- read.csv("./data/multi-omics/DEGs_spatialDRG-general.csv")

head(DEGs)
head(mean_mat)

mat <- mean_mat[rownames(mean_mat) %in% DEGs$gene, ]
mat <- mat[complete.cases(mat), ]
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex),
               point.size = 5,
               ellipse = FALSE,
               ellipse.prob = 0.95,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)
g1 <- g1 + ggtitle("sPCA on spatial-seq DEGs")

print(g1)

pdf(file = paste(PATH_results, "sPCA_spatialDEGs.pdf", sep=""), width = 4, height = 4)
print(g1)
dev.off()

################################################################################

DEGs <- read.csv("./data/multi-omics/DEGs_bulkDRG_fromLGI1analysis.csv")

head(DEGs)
head(mean_mat)

mat <- mean_mat[rownames(mean_mat) %in% DEGs$symbol, ]
mat <- mat[complete.cases(mat), ]
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex),
               point.size = 5,
               ellipse = FALSE,
               ellipse.prob = 0.95,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)
g1 <- g1 + ggtitle("sPCA on bulk-seq DEGs")

print(g1)

pdf(file = paste(PATH_results, "sPCA_bulkDEGs.pdf", sep=""), width = 4, height = 4)
print(g1)
dev.off()

################################################################################

DEGs_a <- read.csv("./data/multi-omics/DEGs_bulkDRG_fromLGI1analysis.csv")
DEGs_b <- read.csv("./data/multi-omics/DEGs_spatialDRG-general.csv")

DEGs <- c(DEGs_a$symbol, DEGs_b$gene)
head(DEGs)

mat <- mean_mat[rownames(mean_mat) %in% DEGs, ]
mat <- mat[complete.cases(mat), ]
pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex),
               point.size = 5,
               ellipse = FALSE,
               ellipse.prob = 0.95,
               labels = NULL,
               labels.size = 4, alpha = 1, var.axes = FALSE,
               circle  = FALSE,
               varname.size = 3,
               varname.adjust = 1.5,
               varname.abbrev = FALSE)

g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)
g1 <- g1 + ggtitle("sPCA on RNA-seq DEGs")

print(g1)

pdf(file = paste(PATH_results, "sPCA_rna-seq-DEGs.pdf", sep=""), width = 4, height = 4)
print(g1)
dev.off()

### extract eigengenes

# extract eigengenes from this to see what drives effect
eigens <- pca$x[,1] #select PC1 frm prcomp output

eigengenes <- pca$rotation[,1]
eigengenes <- eigengenes[order(eigengenes, decreasing=TRUE)]
plot(eigengenes)

# Convert the named vector to a dataframe
eigengenes <- data.frame(
  Gene = names(eigengenes),
  Value = as.numeric(eigengenes)
)

eigengenes <- eigengenes %>%
  arrange(Value) %>%
  mutate(Rank = row_number())

# Identify the top and bottom 5%
n <- nrow(eigengenes)
top_5_percent <- ceiling(n * 0.05)
bottom_5_percent <- top_5_percent

# Add a label column
eigengenes <- eigengenes %>%
  mutate(Label = case_when(
    Rank <= bottom_5_percent ~ Gene,
    Rank > (n - top_5_percent) ~ Gene,
    TRUE ~ NA_character_
  ))

g <- ggplot(eigengenes, aes(x = Rank, y = Value)) +
  geom_point(color = "black") 
g <- g + ggrepel::geom_text_repel(aes(label = Label), 
                                  vjust = -0.5, hjust = 1, max.overlaps = 25, size = 4,
                                  box.padding = 0.35, point.padding = 0.5) 
g <- g + labs(title = "Eigengenes",
              x = "Rank",
              y = "Eigengene Value") 
g <- g + theme_minimal()

print(g)

pdf(paste0(PATH_results, "/ModuleEigengene_DEGs_ranked.pdf"), height = 4, width = 4)
print(g)
dev.off()

# sex comparison
eigengenes <- pca$rotation[,1]
eigengenes <- eigengenes[order(eigengenes, decreasing=TRUE)]

input <- as.data.frame(eigengenes)
input$gene <- rownames(input)
input <- input[c("gene", "eigengenes")]
write.csv(input, file = paste0(PATH_results, "/DEGs_PC1eigens.csv"))

df <- data.frame(eigen=eigens,
                 Sex = colData$Sex)

teigens <- df %>% group_by(Sex, eigen)

bartlett.test(eigen ~ Sex, data = teigens) #test for equal variance (fails, use Welch's, not aov)

df <- as.data.frame(teigens)


oneway.test(eigen ~ Sex, data = df, var.equal = FALSE)

##### change data for plotting (3d v 4w)
g <- ggplot(data=teigens, aes(x=interaction(Sex), y=eigen)) 
g <- g + geom_boxplot(coef = 1.5, aes(fill=Sex)) 
# g <- g + scale_fill_viridis_d() #option = "mako"
g <- g + theme_bw() + theme(
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.text.x = element_text(size=10, angle = 45, hjust= 1),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank())
g <- g + geom_hline(yintercept = 0)

print(g)

pdf(paste0(PATH_results, "/ModuleEigengene_DEGs_boxplot.pdf"), height = 3, width = 3)
print(g)
dev.off()

################################################################################

# # DEGs_a <- read.csv("./data/multi-omics/DEGs_bulkDRG_fromLGI1analysis.csv")
# # DEGs_b <- read.csv("./data/multi-omics/DEGs_spatialDRG-general.csv")
# 
# DEGs <- c(DEGs_a$symbol, DEGs_b$gene, filtered_DARs$gene)
# head(DEGs)
# length(DEGs)
# 
# mat <- mean_mat[rownames(mean_mat) %in% DEGs, ]
# mat <- mat[complete.cases(mat), ]
# pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
# 
# g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
#                groups = interaction(sex),
#                point.size = 5,
#                ellipse = FALSE,
#                ellipse.prob = 0.95,
#                labels = NULL,
#                labels.size = 4, alpha = 1, var.axes = FALSE,
#                circle  = FALSE,
#                varname.size = 3,
#                varname.adjust = 1.5,
#                varname.abbrev = FALSE)
# 
# g1 <- g1 + theme(legend.position = 'right', aspect.ratio=1)
# g1 <- g1 + ggtitle("sPCA on combined multi-omic")
# 
# print(g1)
# 
# pdf(file = paste(PATH_results, "sPCA_on-mulitomic.pdf", sep=""), width = 4, height = 4)
# print(g1)
# dev.off()

#################################################################################

ray_means <- read.csv("../lgi1/output/results_pain_all.csv", header = TRUE)
ray_means <- data.frame(symbol = ray_means$symbol, 
                        rna    = ray_means$baseMean)

mean_mat$protein <- rowMeans(mean_mat)
mean_mat$symbol  <- rownames(mean_mat)
head(mean_mat)

merged_df <- merge(mean_mat, ray_means, by = "symbol")
head(merged_df)

ggplot(merged_df, aes(x = log(rna), y = protein)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a regression line
  labs(title = "Correlation between RNA and protein",
       x = "RNA",
       y = "protein") +
  theme_minimal()

filtered_df <- merged_df[merged_df$rna >= 25, ]

# Calculate the R-squared value
model <- lm(protein ~ log2(rna), data = filtered_df)
r_squared <- summary(model)$r.squared

g <- ggplot(filtered_df, aes(x = log(rna), y = protein)) 
g <- g + geom_point(color = "grey")
g <- g + geom_smooth(method = "lm", se = FALSE, color = "#003391d9") # Add a regression line
g <- g +  labs(title = "Correlation between RNA and protein",
               x = "RNA",
               y = "protein") 
g <- g + annotate("text", x = Inf, y = 10, label = paste("R² = ", round(r_squared, 2)), 
                  hjust = 1.1, vjust = 1.5, size = 5, color = "#003391d9")
g <- g + theme_bw()

print(g)

################################################################################

# check qnTPMs instead
ray_means <- read.csv("../lgi1/output/results_pain_all.csv", header = TRUE)

gene_matching <- data.frame(gene_id = ray_means$X,
                            symbol  = ray_means$symbol)

ray_tpm   <- read.table("../hdrg_immune/processing/ray2023-qnTPMs_subsetforpwoercalcs.txt", header=TRUE)
head(ray_tpm)

ray_tpm <- merge(ray_tpm, gene_matching, by = "gene_id")
ray_tpm$gene_id <- NULL
ray_tpm$rna <- rowMeans(ray_tpm[,1:15])

merged_df <- merge(mean_mat, ray_tpm, by = "symbol")
head(merged_df)

filtered_df <- merged_df[merged_df$rna     >= 10, ]

# Calculate the R-squared value
model <- lm(protein ~ rna, data = filtered_df)
r_squared <- summary(model)$r.squared

# Plot the correlation with R-squared value
g <- ggplot(filtered_df, aes(x = log(rna), y = protein)) 
g <- g + geom_point(color = "grey")
g <- g + geom_smooth(method = "lm", se = FALSE, color = "#003391d9") # Add a regression line
g <- g +  labs(title = "Correlation between RNA and protein",
       x = "RNA",
       y = "protein") 
g <- g + annotate("text", x = Inf, y = 10, label = paste("R² = ", round(r_squared, 2)), 
           hjust = 1.1, vjust = 1.5, size = 5, color = "#003391d9")
g <- g + theme_bw()

print(g)

pdf(file = paste(PATH_results, "correlation-RNAvprotein.pdf", sep=""), width = 5, height = 4)
print(g)
dev.off()

################################################################################

# somalogic correlation
soma  <- read.csv("./output/somalogic/ DEP-analysis-limma-somalogic-healthy.csv", header = TRUE)
input <- read.csv("./output/pilot/ganglia/ DEP-analysis-limma.csv", header = TRUE)

input$dia <- input$logFC
input$EntrezGeneSymbol <- input$X

soma$soma <- soma$logFC

soma  <- soma[abs(soma$logFC) >= 0.2, ]
input <- input[abs(input$logFC) >= 0.2, ]

merged_df <- merge(input, soma, by = "EntrezGeneSymbol")
head(merged_df)



# Calculate the R-squared value
model <- lm(dia ~ soma, data = merged_df)
r_squared <- summary(model)$r.squared

# Plot the correlation with R-squared value
g <- ggplot(merged_df, aes(x = soma, y = dia)) 
g <- g + geom_point(color = "grey")
g <- g + geom_smooth(method = "lm", se = FALSE, color = "#003391d9") # Add a regression line
g <- g +  labs(title = "Sex differences across datasets",
               x = "Somalogic LFC",
               y = "dia-PASEF LFC") 
g <- g + annotate("text", x = Inf, y = Inf, label = paste("R² = ", round(r_squared, 2)), 
                  hjust = 1.1, vjust = 1.5, size = 5, color = "#003391d9")
g <- g + theme_bw()

print(g)

pdf(file = paste(PATH_results, "correlation-somalogicLFC-healthy-filtered.pdf", sep=""), width = 5, height = 4)
print(g)
dev.off()

################################################################################