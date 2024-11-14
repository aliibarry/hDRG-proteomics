library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(viridis)

PATH_results = "./output/"

#-------------------------------------------------------------------------------

pg_types <- readxl::read_excel("./data/multi-omics/interactome_list_v3.1_large.xlsx")
pg_types <- data.frame(Rec_symbol = pg_types$Rec_symbol,
                       Rec_type   = pg_types$Rec_type)

pg_types <- pg_types[!duplicated(pg_types),]

tcounts <- read.csv("./data/multi-omics/ray2023-qnTPMs.csv", row.names = 1)
tcounts$ensembl_gene_id <- rownames(tcounts)

ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_data <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = rownames(tcounts), mart = ensembl)
gene_data <- gene_data[!duplicated(gene_data[,c('ensembl_gene_id')]),]

tcounts <- merge(tcounts, gene_data, by = "ensembl_gene_id")
tcounts <- tcounts %>% distinct(hgnc_symbol, .keep_all = TRUE)
rownames(tcounts) <- tcounts$hgnc_symbol
# tcounts$ensembl_gene_id <- NULL
# tcounts$hgnc_symbol <- NULL

rnaseq_genes <- tcounts$hgnc_symbol
rnaseq_genes <- rnaseq_genes[rnaseq_genes %in% pg_types$Rec_symbol]
rnaseq_genes <- data.frame(genes = rnaseq_genes[!duplicated(rnaseq_genes)])

mat          <- read.csv("./data/processed/matrix-for-limma.csv", header = TRUE, check.names = FALSE)
prot_genes   <- mat[mat$genes %in% pg_types$Rec_symbol, ]

type_rna  <- pg_types[pg_types$Rec_symbol %in% rnaseq_genes$genes, ]
type_prot <- pg_types[pg_types$Rec_symbol %in% mat$genes, ]

#-------------------------------------------------------------------------------

table(type_rna$Rec_type)
table(type_prot$Rec_type)

# Create data frames for RNA and Protein counts
rna_counts <- as.data.frame(table(type_rna$Rec_type))
colnames(rna_counts) <- c("Rec_type", "RNA_Count")

prot_counts <- as.data.frame(table(type_prot$Rec_type))
colnames(prot_counts) <- c("Rec_type", "Protein_Count")
combined_counts  <- merge(rna_counts, prot_counts, by = "Rec_type", all = TRUE)
combined_counts[is.na(combined_counts)] <- 0

combined_counts <- combined_counts %>%
  mutate(Proportion = Protein_Count / RNA_Count)

combined_counts <- combined_counts %>%
  arrange(Proportion) %>%
  mutate(Rec_type = factor(Rec_type, levels = Rec_type))

ggplot(combined_counts, aes(x = Rec_type, y = Proportion, fill = RNA_Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Protein to RNA Counts by Rec_type",
       x = "Rec_type",
       y = "Protein / RNA Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("white", "lightblue", "darkblue"), 
                       trans = "log2", 
                       name = "RNA Count (log scale)")

filtered_counts <- combined_counts %>%
  filter(RNA_Count > 15) %>%
  arrange(Proportion) %>%
  mutate(Rec_type = factor(Rec_type, levels = Rec_type))

ggplot(filtered_counts, aes(x = Rec_type, y = Proportion, fill = RNA_Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Receptor Type Proportions",
       x = "Rec_type",
       y = "Proteomics / Bulk RNA-seq") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_viridis(option = "mako", begin = 0.1, end = 0.9)

pdf(paste0(PATH_results, "receptor_proportions.pdf"), height = 4, width = 6)
ggplot(filtered_counts, aes(x = Rec_type, y = Proportion, fill = RNA_Count)) +
  geom_bar(stat = "identity") +
  labs(title = "Receptor Type Proportions",
       x = "Rec_type",
       y = "Proteomics / Bulk RNA-seq") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_viridis(option = "mako", begin = 0.1, end = 0.9)
dev.off()

#-------------------------------------------------------------------------------

mat <- mat %>%
  distinct(genes, .keep_all = TRUE)

mat_avg_ganglia <- mat %>%
  rowwise() %>%
  mutate(mean_prot = mean(c_across(starts_with("sub-UTD-DN") & ends_with("ganglia")), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(genes, mean_prot)

tcounts_avg <- tcounts %>%
  group_by(hgnc_symbol) %>%
  summarise(mean_rna = mean(c_across(starts_with("X")), na.rm = TRUE)) %>%
  ungroup()

merged_data <- merge(mat_avg_ganglia, tcounts_avg, by.x = "genes", by.y = "hgnc_symbol", all = TRUE)
merged_data <- merge(merged_data, pg_types, by.x = "genes", by.y = "Rec_symbol", all = TRUE)

cor_res <- merged_data %>%
  filter(!is.na(mean_prot) & !is.na(mean_rna)) %>%  # Filter complete cases
  group_by(Rec_type) %>%
  summarise(correlation = cor(mean_prot, mean_rna, use = "complete.obs"), .groups = "drop")

filtered <- combined_counts %>%
  filter(Protein_Count > 15) %>%
  arrange(Proportion) %>%
  mutate(Rec_type = factor(Rec_type, levels = Rec_type))

cor_res <- as.data.frame(cor_res[cor_res$Rec_type %in% filtered$Rec_type, ])

rownames(cor_res) <- cor_res$Rec_type
cor_res$Rec_type  <- NULL
cor_mat <- as.matrix(cor_res)

head(cor_res)

pdf(paste0(PATH_results, "receptor_correlations.pdf"), height = 5, width = 4)
Heatmap(cor_mat,
        col=viridis(100),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_title = "",
        row_dend_side = "left")
dev.off()
