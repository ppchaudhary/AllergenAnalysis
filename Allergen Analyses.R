library(tidyverse)
library(dplyr)
library(pheatmap)
library(ggplot2)
#install.packages("readxl")
library(readxl)
#install.packages("ggpubr")
library(ggpubr)

setwd("Desktop/ALEX2/")

#Reading In and Cleaning Data----
reports <- list.files(path = "AssayReports", pattern = "\\.csv$", full.names = TRUE) #compiling file names as a list

report_list <- list() #creating list for file names to go in from loop

for (report in reports) { #loop 
  data <- read.csv(report, sep = ";") #reading in files and specifying delimiter
  report_list[[report]] <- data } #adding read in files to list

mergeddf <- bind_rows(report_list)

mergeddf <- mergeddf[!(mergeddf$AllergenName %in% "GD1"),] #removing the guide dots from dataframe
mergeddf <- mergeddf[!(mergeddf$AllergenName %in% "GD2"),]
mergeddf <- mergeddf[!(mergeddf$AllergenName %in% "GD3"),]
mergeddf <- mergeddf[!(mergeddf$AllergenName %in% "GD4"),]

mergeddf <- mergeddf[!duplicated(mergeddf), ] #removing rows that were duplicated 

#Creating a File with Names for Each Allergens
allergynames <- read_excel("export_2025-02-24_16-21-22.xlsx", range = "A15:B315")
colnames(allergynames) <- c("Allergen", "Name")

allergynames <- allergynames %>% 
  add_row(Allergen = "tIgE", Name = "Total IgE") #creating a name for tIgE

mergenames <- merge(mergeddf, allergynames, by.x = "AllergenName", by.y = "Allergen", all.x = TRUE)

mergenames <- mergenames %>%
  mutate(Group = substr(mergenames$SampleCode, 1L, 2L)) %>% #creates group column based on AD, TSW, or HV
  select("Group", everything())

mergenames[mergenames == 'PI'] <- 'TSW'
mergenames <- mergenames %>%
  mutate_all(funs(str_replace(., "< 0.10", "0"))) %>%
  mutate_all(funs(str_replace(., "< 20.00", "20.00"))) %>%
  mutate_all(funs(str_replace(., "> 2500.00", "2500.00")))

df <- mergenames %>% 
  select(Name, AllergenName, Group, CalibratedValue) %>% #this selects only the specified rows
  group_by(Name)

sum(is.na(df)) #checking how many NA's are in df because errors were popping up

which(is.na(df), arr.ind=TRUE) #finding which rows have NA, saw it was tIgE 
#and GD rows so removed GD and added name for tIgE

df[4] <- lapply(df[4], as.numeric) #converting "CalibratedValue" col to numeric 
#because errors while determining mean

class(df$CalibratedValue) #making sure CalibratedValue column is numeric

#summarising df when grouped by Allergen Name and specific epitopes
test1 <- df %>% 
  group_by(Name, Group) %>%
  summarise(
    count = sum(CalibratedValue != 0),
    mean = mean(CalibratedValue, na.rm = TRUE),
    sd = sd(CalibratedValue, na.rm = TRUE))

test2 <- df %>%
  group_by(AllergenName, Group) %>%
  summarise(
    count = sum(CalibratedValue != 0),
    mean = mean(CalibratedValue, na.rm = TRUE),
    sd = sd(CalibratedValue, na.rm = TRUE))

#Lines 78-101 were checking individual allergens, but created a loop so disregard
#aa
#alternaria <- subset(df, Name == "Alternaria alternata")

#group_by(alternaria, Name, Group) %>%
  #summarise(
    #count = sum(CalibratedValue != 0),
    #mean = mean(CalibratedValue, na.rm = TRUE),
    #sd = sd(CalibratedValue, na.rm = TRUE))

#res.aov <- aov(CalibratedValue ~ Group, data = alternaria)
#summary(res.aov)

#dog
#dog <- subset(df, Name == "Dog")

#group_by(dog, Name, Group) %>%
  #summarise(
    #count = sum(CalibratedValue != 0),
    #mean = mean(CalibratedValue, na.rm = TRUE),
    #sd = sd(CalibratedValue, na.rm = TRUE))

#res.aov1 <- aov(CalibratedValue ~ Group, data = dog)
#summary(res.aov1) #significant difference in dog allergy
#TukeyHSD(res.aov1)

#For Loop with One-Way Anova with Groups----
allergens <- test1[["Name"]] #creating a list with all of the values form the "Name" column
allergens <- allergens[!duplicated(allergens)] #removing duplicates

allergen_list <- list()

for (allergen in allergens) { #loop
  x <- subset(df, Name == allergen) 
  allergen_list[[allergen]] <- x  #this loop is creating multiple dataframes for each allergen
group_by(x, Name, Group) %>%
  summarise( #this is computing summary statistics for each one
    count = n(),
    mean = mean(CalibratedValue, na.rm = TRUE),
    sd = sd(CalibratedValue, na.rm = TRUE)
  )
y <- aov(CalibratedValue ~ Group, data = x) #performing the anova
result <- summary(y) #storing the anova as a variable
print(paste("Summary of One Way Anova for", allergen))
print(result)
}

#this if statement wasn't working
#if ((results)[[1]]["Pr(>F)"] <= 0.05) {
#print(paste("Pairwise comparisons with Tukey test:"))
#print(TukeyHSD(model, conf.level = 0.95))}

#For Loop with One-Way Anova with Specific Allergen
names <- test2[["AllergenName"]] #creating a list with all of the values form the "Name" column
names <- names[!duplicated(names)] #removing duplicates

name_list <- list()

for (name in names) { #loop
  x <- subset(df, AllergenName == name) 
  name_list[[name]] <- x  #this loop is creating multiple dataframes for each allergen
  group_by(x, AllergenName, Group) %>%
    summarise( #this is computing summary statistics for each one
      count = sum(CalibratedValue !=0),
      mean = mean(CalibratedValue, na.rm = TRUE),
      sd = sd(CalibratedValue, na.rm = TRUE)
    )
  y <- aov(CalibratedValue ~ Group, data = x) #performing the anova
  results <- summary(y) #storing the anova as a variable
  print(paste("Summary of One Way Anova for", name))
  print(results)
}

#Source:https://www.sthda.com/english/wiki/one-way-anova-test-in-r
plot <- ggboxplot(dog, x = "Group", y = "CalibratedValue", 
          color = "Group", palette = c("red", "blue", "orange"),
          order = c("HV", "AD", "TSW"),
          ylab = "IgE Reactivity to Dog Allergens")
plot


plot1 <- ggline(dog, x = "Group", y = "CalibratedValue", 
                add = c("mean_se", "jitter"), 
                order = c("HV", "AD", "TSW"),
                ylab = "IgE Reactivity to Dog Allergens", xlab = "Group")

plot1

plot(res.aov1, 1)

#need to write code to make plots for all significant results

#does it make sense to test allergy with multiple allergens (Dog) vs multiple allergens (Can f_Fd1)?
#Ian said yes for aeroallergens, but not for food