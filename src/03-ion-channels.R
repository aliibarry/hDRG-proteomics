library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(msigdbr)

PATH_results = "./output/"
mat     <- read.csv("./data/processed/matrix-for-limma.csv", header = TRUE, check.names = FALSE)
colData <- read.csv("./data/processed/colData-for-limma.csv")

################################################################################

mat <- mat %>% distinct(genes, .keep_all = TRUE)

goi <- read.table("./data/ionchannels.csv", fill = TRUE, header = TRUE, sep=",")
goi <- goi[which(goi$HGNC.symbol %in% mat$genes==TRUE), ]

df <- mat %>% pivot_longer(cols = -c(proteins, genes),
               names_to = "mergedID",   
               values_to = "expression") 

df <- df %>% left_join(colData, by = "mergedID") 

df <- df[df$genes %in% goi$HGNC.symbol, ]

df_means <- df %>%
  group_by(genes, Tissue, Sex) %>%
  summarise(expression = mean(expression, na.rm = TRUE), .groups = 'drop')

# select ion channels
df_channels <- df_means %>% 
  filter(str_detect(genes, "^TRP") | 
          str_detect(genes, "^KCN") | 
          str_detect(genes, "^SCN") |
          str_detect(genes, "^GAB") | 
          str_detect(genes, "^ASIC") | 
          str_detect(genes, "^CNG") | 
          str_detect(genes, "^RYR") | 
          str_detect(genes, "^CAC") | 
          str_detect(genes, "^HVC"))

df_gpcr <- df_means %>% 
  filter(
    str_detect(genes, "^ADRA") | 
      str_detect(genes, "^ADRB") | 
      str_detect(genes, "^OPR") | 
      str_detect(genes, "^GPR") | 
      str_detect(genes, "^SSTR") | 
      str_detect(genes, "^DRD") | 
      str_detect(genes, "^HTR") | 
      str_detect(genes, "^CNR") | 
      str_detect(genes, "^EDNR") | 
      str_detect(genes, "^FPR") | 
      str_detect(genes, "^MC") | 
      str_detect(genes, "^TAS") | 
      str_detect(genes, "^AVPR") | 
      str_detect(genes, "^P2RY") | 
      str_detect(genes, "^LGR")
  )


# Define a function to create plots
plot_genes <- function(data, gene_pattern) {
  ggplot(data, aes(x = genes, y = interaction(Sex))) +
    facet_grid(Tissue~.) +
    scale_colour_viridis_c(option = "mako", end = .90, limits = c(8, 18)) + 
    geom_point(aes(col = expression, size = expression)) +
    scale_size_continuous(limits = c(8, 18)) +  # Set limits for the point size
    theme_bw() +
    labs(title = paste(gene_pattern, "proteins")) +
    theme(panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )
  }

# Create and display the plots
plot_channels <- plot_genes(df_channels, "Ion channel")
plot_gpcr     <- plot_genes(df_gpcr, "GPCR")

plot_channels
plot_gpcr

pdf(paste0(PATH_results, "Channels.pdf"), height = 3, width = 7)
plot_channels
dev.off()

pdf(paste0(PATH_results, "GPCR.pdf"), height = 3, width = 5)
plot_gpcr
dev.off()

################################################################################

# calculate proportion of proteins detected in each pain-related gene set

mat <- read.csv("./data/processed/matrix-for-limma.csv", header = TRUE, check.names = FALSE)
mat <- mat %>% distinct(genes, .keep_all = TRUE)

pain_genes   = read.csv("../pain-classifier/data/pg.csv", row.names = 1)
immune_genes = read.csv("../hdrg_proteomics/data/immunedisease_gs.csv")
GO_gene_sets = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

################################################################################

gene.ids <- data.frame(symbol = mat$genes)
head(gene.ids)

gs_of_interest <- c("HSNPS", "PL", "PGD",
                    "sensory_pain", "response_to_pain")

pain_genes <- pain_genes[pain_genes$gs %in% gs_of_interest,]

# #other options
# gs_of_interest <- c("pain_drugs", "chronic_pain_drugs", "neuropathic_pain_drugs")
# pain_genes <- immune_genes

bonus <- GO_gene_sets[GO_gene_sets$gs_name == "GOBP_INFLAMMATORY_RESPONSE", ]
bonus <- bonus[, c("human_ensembl_gene", "gs_name", "gene_symbol")]
bonus$gs <- "Inflammatory response"
bonus$gene <- bonus$human_ensembl_gene
bonus$symbol <- bonus$gene_symbol

bonus <- as.data.frame(bonus)
bonus <- bonus[, c("symbol", "gene", "gs")]

#to add inflammatory_response
pain_genes <- rbind(pain_genes, bonus)

gene.ids$proteomics <- 1
pain_genes$gene_set <- 1
pain_genes$gene <- NULL

head(gene.ids)
#gene.ids <- gene.ids %>% dplyr::distinct()

proportions <- merge(gene.ids, pain_genes, by = "symbol", all.y = TRUE)
proportions[is.na(proportions)] <- 0
head(proportions)

proportions <- proportions %>%
  dplyr::group_by(gs) %>%
  dplyr::summarise(
    prop = mean(as.integer(proteomics)),
    count = sum(as.integer(proteomics == 1))
  )

head(proportions)

g <- ggplot(proportions, aes(x = gs, y = prop, fill = log(count))) 
g <- g + geom_bar(stat = "identity") + 
  scale_fill_gradient(low = "#e8740036", high = "#003391d9",
                      breaks = c(20, 100, 200),
                      labels = c(20, 100, 200))
g <- g + labs(y = "Proportion", size = 15)
g <- g + theme_bw() + ggtitle("Pain Pathway Overlap") +
  theme(axis.text.y = element_text(size= 12, colour= "black", hjust = 1), 
        axis.text.x = element_text(size=12, angle = 45, hjust= 1), 
        legend.text = element_text(size=12), 
        axis.title.x = element_blank(),
        plot.title=element_text(size=rel(1), hjust = 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g)

pdf(paste0("boxplot_immune.pdf"), height = 4, width = 3)
print(g)
dev.off()


