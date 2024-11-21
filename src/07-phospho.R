library(dplyr)
library(ggplot2)
library(readxl)

dir.create("./output/phospho/")
PATH_results = "./output/phospho/"

array_df <- readxl::read_excel("./data/processed/phosphoarray_sourcedata.xlsx")

arraya <- array_df[array_df$Array %in% "A",]

tail(arraya)

pbs_values <- arraya %>%
  filter(Pair == "PBS") %>%
  group_by(Donor, Round) %>%
  summarize(PBS_Mean = mean(Mean, na.rm = TRUE), .groups = "drop")

adjusted_arraya_summary <- arraya %>%
  filter(Pair != "PBS") %>%
  group_by(Donor, Round, Pair, Sex) %>%
  summarize(Mean_ROI = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  left_join(pbs_values, by = c("Donor", "Round")) %>%
  mutate(Adjusted_Mean_ROI = PBS_Mean - Mean_ROI)

# Calculate average Adjusted_Mean_ROI by Pair and Sex for the bar plot
A_pair_avg_summary <- adjusted_arraya_summary %>%
  group_by(Pair, Sex) %>%
  summarize(Average_Adjusted_Mean_ROI = mean(Adjusted_Mean_ROI, na.rm = TRUE), .groups = "drop")

g <- ggplot(A_pair_avg_summary, aes(x = Pair, y = Average_Adjusted_Mean_ROI, fill = Sex)) 
g <- g + geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) 
g <- g + geom_jitter(data = adjusted_arraya_summary, 
              aes(x = Pair, y = Adjusted_Mean_ROI, color = Sex), 
              width = 0.15, size = 2, shape = 21, fill = "white") 
g <- g + labs(title = "",
       x = "Site",
       y = "Normalized ROI Intensity") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g)

# pdf(paste0(PATH_results, "array-A.pdf"), height = 4, width = 6)
# print(g)
# dev.off()

#-------------------------------------------------------------------------------

arrayb <- array_df[array_df$Array %in% "B",]

head(arrayb)

pbs_values <- arrayb %>%
  filter(Pair == "PBS") %>%
  group_by(Donor, Round) %>%
  summarize(PBS_Mean = mean(Mean, na.rm = TRUE), .groups = "drop")

adjusted_arrayb_summary <- arrayb %>%
  filter(Pair != "PBS") %>%
  group_by(Donor, Pair, Sex, Round) %>%
  summarize(Mean_ROI = mean(Mean, na.rm = TRUE), .groups = "drop") %>%
  left_join(pbs_values, by = c("Donor", "Round")) %>%
  mutate(Adjusted_Mean_ROI = PBS_Mean - Mean_ROI)

# Calculate average Adjusted_Mean_ROI by Pair and Sex for the bar plot
B_pair_avg_summary <- adjusted_arrayb_summary %>%
  group_by(Pair, Round, Sex) %>%
  summarize(Average_Adjusted_Mean_ROI = mean(Adjusted_Mean_ROI, na.rm = TRUE), .groups = "drop")

g <- ggplot(B_pair_avg_summary, aes(x = Pair, y = Average_Adjusted_Mean_ROI, fill = Sex)) 
g <- g + geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) 
g <- g + geom_jitter(data = adjusted_arrayb_summary, 
                     aes(x = Pair, y = Adjusted_Mean_ROI, color = Sex), 
                     width = 0.15, size = 2, shape = 21, fill = "white") 
g <- g + labs(title = "",
              x = "Site",
              y = "Normalized ROI Intensity") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g)

# pdf(paste0(PATH_results, "array-B.pdf"), height = 4, width = 6)
# print(g)
# dev.off()

#-------------------------------------------------------------------------------

mean_summary    <- bind_rows(A_pair_avg_summary, B_pair_avg_summary)
general_summary <- bind_rows(adjusted_arraya_summary, adjusted_arrayb_summary)

g <- ggplot(mean_summary, aes(x = Pair, y = Average_Adjusted_Mean_ROI, fill = Sex)) 
g <- g + geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.7, width = 0.7) 
g <- g + geom_jitter(data = general_summary , 
                     aes(x = Pair, y = Adjusted_Mean_ROI, color = Sex), 
                     width = 0.15, size = 2, shape = 21, fill = "white") 
g <- g + labs(title = "",
              x = "Site",
              y = "Normalized ROI Intensity") +
  #scale_y_log10() +
  scale_y_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::label_number(accuracy = 1)) +
  theme_bw() +
  scale_fill_manual(values = c("F" = "#d9a2e4ff", "M" = "#5abad2ff")) + 
  scale_color_manual(values = c("F" = "#a939bdff", "M" = "#4fa3b7ff")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(g)

# pdf(paste0(PATH_results, "array-combined.pdf"), height = 4, width = 9)
# print(g)
# dev.off()

#-------------------------------------------------------------------------------

g <- ggplot(general_summary, aes(x = Pair, y = Adjusted_Mean_ROI, fill = Sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplot
  geom_jitter(aes(color = Sex), width = 0.2, size = 1.5, alpha = 0.8) +  # Jitter for individual points
  labs(title = " ",
       x = "Pair",
       y = "Normalized ROI Intensity") +
  scale_y_continuous(trans = "log10", 
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::label_number(accuracy = 1)) +
  theme_bw() +
  scale_fill_manual(values = c("F" = "#d9a2e4ff", "M" = "#5abad2ff")) +
  scale_color_manual(values = c("F" = "#a939bdff", "M" = "#4fa3b7ff")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

print(g)

pdf(paste0(PATH_results, "array-boxplot.pdf"), height = 5, width = 12)
print(g)
dev.off()

#-------------------------------------------------------------------------------

library(multcomp)
library(emmeans)
library(dplyr)

general_summary$Sex  <- as.factor(general_summary$Sex)
general_summary$Pair <- as.factor(general_summary$Pair)

# Perform ANOVA
anova_model <- aov(Adjusted_Mean_ROI ~ Sex*Pair, data = general_summary)
summary(anova_model)

posthoc <- TukeyHSD(anova_model)

results <- posthoc$`Sex:Pair`
comparison_names <- rownames(results)

# output pairs of interest
is_within_pair <- function(comparison) {
  levels <- strsplit(comparison, "-")[[1]]
  pair1 <- sub("^[MF]:", "", levels[1])
  pair2 <- sub("^[MF]:", "", levels[2])
  
  # Check if the pairs are identical
  return(pair1 == pair2)
}

# Apply the function to filter within-pair comparisons
within_pair_indices <- sapply(comparison_names, is_within_pair)
results <- results[within_pair_indices, ]

print(results)

#-------------------------------------------------------------------------------

# pair specific posthoc testing

emmeans_results <- emmeans(anova_model, ~ Sex | Pair)

# Run pairwise comparisons within pairs only
within_pairs <- pairs(emmeans_results)
posthoc_df   <- as.data.frame(within_pairs)

posthoc <- subset(posthoc_df, Pair %in% unique(general_summary$Pair))
print(posthoc)
