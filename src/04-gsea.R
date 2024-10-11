library("clusterProfiler")
library(msigdbr)
library(ggplot2)
library(dplyr)

dir.create("./output/tissue/GSEA/")
PATH_results = "./output/tissue/GSEA/"

input <- read.csv("./output/tissue/ DEP-limma-tissue.csv", header = TRUE)
input <- input[c("X", "adj.P.Val", "logFC")]

################################################################################

# load relevant gene sets
all_gene_sets = msigdbr(species = "Homo sapiens")
hallmark_sets = msigdbr(species = "Homo sapiens", category = "H")
GO_gene_sets  = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
#MF_gene_sets  = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
#C8_gene_sets  = msigdbr(species = "Homo sapiens", category = "C8") #cell type markers
#reactome_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")

################################################################################

filtered_dge <- input %>%
  dplyr::arrange(dplyr::desc(abs(logFC))) %>%
  dplyr::distinct(X, .keep_all = TRUE)

lfc_vector <- filtered_dge$logFC
names(lfc_vector) <- filtered_dge$X

# Sort log2 fold change values in descending order
lfc_vector <- sort(lfc_vector, decreasing = TRUE)

set.seed(52)

gsea_results <- GSEA(
  geneList = lfc_vector,
  minGSSize = 25, 
  #maxGSSize = 500,
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(hallmark_sets, gs_name, gene_symbol)
)

gsea_result_df <- data.frame(gsea_results@result)
write.csv(gsea_result_df, file = paste0(PATH_results, "GSEA_results_tissuetypes.csv"))

g <- ggplot(gsea_result_df, aes(x=(Description), y=NES, colour=p.adjust, size=setSize))
g <- g + geom_point() + theme_bw() + ggtitle("hDRG proteomics, tissue type") +
  theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
        axis.text.x = element_text(size=10, angle = 45, hjust= 1), 
        legend.text = element_text(size=10), 
        axis.title.x = element_blank(),
        plot.title=element_text(size=rel(1), hjust = 1)) +
  theme(plot.margin=unit(c(0.3,0.3,0.3,0.3),"cm")) +  
  labs(y="Enrichment Score", colour="p value", size="Count")

pdf(paste0(PATH_results, "GSEA_dots_hallmark_tissuetypes.pdf"), height = 5, width = 4)
print(g)
dev.off()
