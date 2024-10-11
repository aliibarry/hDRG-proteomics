library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)

# use ganglia/NR specific matrices
PATH_results = "./output/"
mat     <- read.csv("./data/processed/matrix-for-limma.csv", header = TRUE, check.names = FALSE) #example on full matrix
colData <- read.csv("./data/processed/colData-for-limma.csv")

################################################################################

# sPCA based on DEGs and equiv for sex in hDRG
DARs <- read.csv("./data/multi-omics/DARs_Franco-Enzastiga2024.csv")
filtered_DARs <- DARs %>% filter(!(Chr %in% c("chrX", "chrY")))

head(filtered_DARs)
head(mat)

mat <- mat %>% distinct(genes, .keep_all = TRUE)
rownames(mat) <-mat$genes

mat$proteins <- NULL
mat$genes    <- NULL

# reorder colData to match matrix
colData <- colData[colData$mergedID %in% colnames(mat), ]

index   <- match(colnames(mat), colData$mergedID)
colData <- colData[index, ]
sex     <- as.factor(colData$Sex) #sex factor as needed
tissue  <- as.factor(colData$Tissue)

mat <- mat[rownames(mat) %in% filtered_DARs$gene, ] #filtered_
mat <- mat[complete.cases(mat), ]

pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

g1 <- ggbiplot(pca, choices = c(1,2), scale = 0,
               groups = interaction(sex, tissue),
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

################################################################################

# Plotting and stats

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

write.csv(input, file = paste0(PATH_results, "/DARs_PC1eigens.csv"))

df <- data.frame(eigen=eigens,
                 Sex = colData$Sex)

teigens <- df %>% group_by(Sex, eigen)

bartlett.test(eigen ~ Sex, data = teigens) #test for equal variance (fails, use Welch's, not aov)

df <- as.data.frame(teigens)

oneway.test(eigen ~ Sex, data = df, var.equal = TRUE)

# BH <- pairwise.t.test(teigens4W$eigen, teigens4W$Pop_Cond,
#                 p.adjust.method = "BH", pool.sd = FALSE)
# 
# pairwise.wilcox.test(teigens$eigen, teigens$Sex,
#                      p.adjust.method = "BH")
# 
# kruskal.test(eigen ~ Sex, data = teigens)

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
