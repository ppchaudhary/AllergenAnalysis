AllergenAnalysis is an R-based pipeline for processing and analyzing IgE reactivity to allergens. It reads assay reports, cleans and integrates data, performs ANOVA to identify significant allergens, and visualizes the results using boxplots and beeswarm plots.

Features
📊 Data Processing: Reads and merges allergen assay reports.
🔍 Filtering & Cleaning: Removes unwanted allergens and duplicates.
📈 Statistical Analysis: Uses ANOVA to detect significant allergens.
🎨 Visualization: Generates boxplots and beeswarm plots for data interpretation.
📁 Automated Output: Saves plots for significant allergens for easy review.
Installation
Prerequisites
Ensure you have R and the required libraries installed:

install.packages(c("tidyverse", "dplyr", "pheatmap", "ggplot2", "readxl", "ggpubr", "ggbeeswarm"))
Usage
Clone the repository
Example data is alsp provided along with the output files.

git clone https://github.com/luckypgi/AllergenAnalysis.git
cd AllergenAnalysis
Run the script in R

source("Edited script.R")
View generated plots

Boxplots and beeswarm plots will be saved in the working directory.
Input & Output
Input: CSV assay reports and an Excel file with allergen names.
Output: Filtered dataset, significant allergen list, and statistical plots.
Contributing
Contributions are welcome! Feel free to submit issues or pull requests.

License
This project is licensed under the MIT License.
