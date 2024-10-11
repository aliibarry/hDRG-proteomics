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

