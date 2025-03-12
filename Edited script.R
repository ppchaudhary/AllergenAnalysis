library(tidyverse)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(readxl)
library(ggpubr)

setwd("Desktop/ALEX2/")

# Read all assay reports
reports <- list.files(path = "AssayReports", pattern = "\\.csv$", full.names = TRUE)
report_list <- lapply(reports, function(report) read.csv(report, sep = ";"))
mergeddf <- bind_rows(report_list)

# Remove unwanted allergens and duplicates
mergeddf <- mergeddf %>% filter(!(AllergenName %in% c("GD1", "GD2", "GD3", "GD4")))
mergeddf <- mergeddf %>% distinct()

# Read allergen names
allergynames <- read_excel("export_2025-02-24_16-21-22.xlsx", range = "A15:B315")
colnames(allergynames) <- c("Allergen", "Name")
allergynames <- allergynames %>% add_row(Allergen = "tIgE", Name = "Total IgE")

# Merge allergen names with the dataset
mergenames <- merge(mergeddf, allergynames, by.x = "AllergenName", by.y = "Allergen", all.x = TRUE)

# Add Group column
mergenames <- mergenames %>%
  mutate(Group = substr(SampleCode, 1L, 2L)) %>%
  select("Group", everything())

# Replace PI with TSW
mergenames[mergenames == 'PI'] <- 'TSW'

# Standardize values
mergenames <- mergenames %>%
  mutate_all(~str_replace(., "< 0.10", "0")) %>%
  mutate_all(~str_replace(., "< 20.00", "20.00")) %>%
  mutate_all(~str_replace(., "> 2500.00", "2500.00"))

# Convert CalibratedValue to numeric
df <- mergenames %>%
  select(Name, AllergenName, Group, CalibratedValue) %>%
  mutate(CalibratedValue = as.numeric(CalibratedValue)) %>%
  group_by(Name)

# Perform ANOVA and filter significant allergens
significant_allergens <- df %>%
  group_split() %>%
  map_df(~ {
    allergen_name <- unique(.x$Name)
    if (n_distinct(.x$Group) > 1) {
      model <- aov(CalibratedValue ~ Group, data = .x)
      p_value <- summary(model)[[1]]$`Pr(>F)`[1]
      if (!is.na(p_value) && p_value < 0.05) {
        tibble(Name = allergen_name, p_value = p_value)
      } else {
        NULL
      }
    }
  })

# Generate and save boxplots for significant allergens
output_dir <- getwd()  # Current working directory

for (i in 1:nrow(significant_allergens)) {
  allergen <- significant_allergens$Name[i]
  p_val <- significant_allergens$p_value[i]
  
  x <- df %>% filter(Name == allergen)
  
  plot <- ggboxplot(x, x = "Group", y = "CalibratedValue",
                    color = "Group", palette = c("red", "blue", "orange"),
                    order = c("HV", "AD", "TSW"),
                    ylab = paste("IgE Reactivity to", allergen),
                    title = paste("Boxplot for", allergen)) +
    annotate("text", x = 1.5, y = max(x$CalibratedValue, na.rm = TRUE),
             label = paste("ANOVA p-value:", format(p_val, digits = 3)),
             size = 5, color = "black")
  
  # Save the plot
  ggsave(filename = paste0(output_dir, "/", allergen, "_significant_plot.png"),
         plot = plot, width = 6, height = 4)
}




####ggbeeswarm plots
library(ggbeeswarm)

ggplot(df %>% filter(Name %in% significant_allergens$Name),
       aes(x = Group, y = CalibratedValue, color = Group)) +
  geom_beeswarm(size = 2, dodge.width = 0.8) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  facet_wrap(~ Name, scales = "free") +
  labs(title = "Beeswarm Plot of Significant Allergens",
       y = "IgE Reactivity") +
  theme_minimal()



