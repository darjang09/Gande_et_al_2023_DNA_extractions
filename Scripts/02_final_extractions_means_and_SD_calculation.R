## Calculations for additional information regarding the data (Table S1 and S2)

library(data.table)
library(ggplot2)
library(outliers)
library(car)
library(rcompanion)
library(FSA)
library(dplyr)


## Table S2

# Import the data
df <- read.table("../../All_the_data/06_Results/01_extraction_optimization/Extractions_2022/Statistics/data_for_stats_minerals_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

# calculate means

means <- tapply(df$Value, list(df$Mineral, df$Buffer, df$Extraction), mean)

means_table <- as.data.frame(round(means, 2))

# calculate SD

stdev <- tapply(df$Value, list(df$Mineral, df$Buffer, df$Extraction), sd)

sd_table <- as.data.frame(round(stdev, 2))


### Table S1 - unbound DNA

df <- read.table("../../All_the_data/06_Results/01_extraction_optimization/Extractions_2022/Statistics/data_for_stats_unbound_minerals_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

# calculate means

means <- tapply(df$Value, list(df$Mineral, df$Buffer, df$Extraction), mean)

means_table <- as.data.frame(round(means, 2))

# calculate SD

stdev <- tapply(df$Value, list(df$Mineral, df$Buffer, df$Extraction), sd)

sd_table <- as.data.frame(round(stdev, 2))
