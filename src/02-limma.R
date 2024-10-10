library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("./output/tissue")
PATH_results="./output/tissue/"

expression <- read.csv("./data/processed/matrix-for-limma.csv", row.names = 1, check.names = FALSE)
colData    <- read.csv("./data/processed/colData-for-limma.csv")

################################################################################

rownames(colData) <- colData$mergedID

# test unique genes
expression <- expression %>%
  distinct(genes, .keep_all = TRUE)

row.names(expression) <- expression$genes
expression$genes <- NULL

# filter for 80% present
NAs <- rowMeans(is.na(expression))
expression <- expression[NAs <= 0.8, ]
dim(expression)

tissue <- as.factor(colData$Tissue)
design <- model.matrix(~ 0 + Tissue, data = colData)

head(design)

fit <- lmFit(expression, design)

contrast.matrix <- makeContrasts(nerve.root_vs_ganglia = Tissuenerve.root - Tissueganglia, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust="BH", number=Inf)
head(results)

sig_de <- results[results$adj.P.Val < 0.05, , drop = FALSE]
sig_de <- sig_de[abs(sig_de$logFC) > 1, ]

head(sig_de)

write.csv(results, paste(PATH_results, "DEP-limma-tissue.csv"))
write.csv(sig_de, paste(PATH_results,  "DEP-limma-tissue_significant.csv"))

################################################################################

expression$genes <- row.names(expression)

df <- expression %>%
  pivot_longer(cols = -c(genes),
               names_to = "mergedID",   
               values_to = "expression") 

df      <- df %>% left_join(colData, by = "mergedID") 
df.filt <- df %>% filter(genes %in% c("TRPV1", "MBP"))

g <- ggplot(df.filt, aes(x = Tissue, y = expression, fill = Tissue)) 
g <- g + geom_boxplot() 
g <- g +  facet_wrap(~genes) + theme_bw()
g <- g +labs(title = "", y = "log2(expression)", x = "Sex") 
g <- g + scale_fill_manual(values = c("F" = "#914c83", "M" = "grey"))
g <- g + scale_fill_manual(values = c("ganglia" = "#914c83", "nerve.root" = "grey"))

pdf(file = paste(PATH_results, "DEP_example.pdf", sep=""), width = 5, height = 3)
print(g)
dev.off()

################################################################################

# volcano plotting, modified from Georgios Baskozos, University of Oxford
res       <- as.data.frame(results) 
mutateddf <- mutate(res, Sig=ifelse(res$adj.P.Val<0.05 & abs(res$logFC)>1, "FDR < 0.05 & LFC > 1", ifelse("NS")))
input     <- cbind(gene=rownames(res), mutateddf) 

volc = ggplot(input, aes(logFC, -log10(P.Value))) + geom_point(aes(col=Sig)) +
  #scale_color_manual(values = c("#0D0887FF","#9512A1FF", "grey")) + 
  scale_color_manual(values = c("#B63679ff", "grey")) + 
  ggrepel::geom_text_repel(data=subset(input, input$gene %in% rownames(sig_de)),
                           aes(label=gene), size=4, segment.alpha= 0.2, force =0.5, max.overlaps=20) 
volc <- volc + theme_bw() + theme(aspect.ratio=1)
volc <- volc + theme(legend.position="bottom", axis.text.y = element_text(size= 12, face="bold"), 
                     axis.title.y = element_text(size=14), axis.title.x = element_text(size= 14), 
                     axis.text.x = element_text(size= 12), legend.title=element_text(size=14), 
                     legend.text=element_text(size=14), plot.title=element_text(size=12, hjust = 0.5)) 
#+ ggtitle("hDRG Proteomics")

pdf(file = paste0(PATH_results, "volcano.pdf"), height = 4, width = 4)
print(volc)
dev.off()

################################################################################

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)

# library("msigdbr")
# msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
#   dplyr::select(gs_name, entrez_gene) %>%
#   dplyr::rename(ont = gs_name, gene = entrez_gene)
# 
# head(msig_h)
# 
# gene.df <- data.frame(gene = rownames(results))
# 
# gene.df <- bitr(gene.df$gene, fromType = "SYMBOL",
#                 toType = c("ENSEMBL", "SYMBOL"),
#                 OrgDb = org.Hs.eg.db)
# 
# DEPs <- gene.df[gene.df$SYMBOL %in% rownames(sig_de), ]

upregulated   <- rownames(sig_de[sig_de$logFC > 1, ])
downregulated <- rownames(sig_de[sig_de$logFC < 1, ])

ego <- enrichGO(gene          = downregulated, #swap as needed
                universe      = rownames(results),
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
goplot(ego)


pdf(file = paste(PATH_results, "GO-tissue_nerve-barplot.pdf", sep=""), width = 5, height = 3)
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")
dev.off()

pdf(file = paste(PATH_results, "GO-tissue_nerve-upset.pdf", sep=""), width = 7, height = 4)
upsetplot(ego)
dev.off()

write.csv(ego, paste(PATH_results, "GO-tissue_nerve.csv"))
