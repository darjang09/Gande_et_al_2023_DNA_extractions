## Script for plotting Supplementary Figure S2 - testing isolation methods without mineral component

## Load packages

library(data.table)
library(ggplot2)
library(outliers)
library(car)
library(rcompanion)
library(FSA)
library(dplyr)
library(scales)
library(plyr)
library(patchwork)

# Import the data
df <- read.table("data_for_stats_controls_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

## Process the data
agg=aggregate(Value~Buffer*Extraction, data=df, FUN="mean") #mean
agg$sd=aggregate(Value~Buffer*Extraction, data=df, FUN="sd")$Value #add the SD
agg$se <- with(agg, sd/2) #add the standard error (SD/2)

df.plot <- agg

names(df.plot)[c(1:4)] <- c("Buffer","Extraction","Value","SD")
str(df.plot)
df.plot$Buffer <- factor(df.plot$Buffer, levels=c("Input_DNA", "Mag_Beads_only", "Amicon + Mag_Beads", "24h in ASW + Amicon + Mag_Beads",
                                                      "Input_DNA", "Precipitation", "PCI + Precipitation", "Amicon + Precipitation", "PCI + Amicon + Precipitation"),
                           labels=c("Input_DNA", "Mag_Beads_only", "Amicon + Mag_Beads", "24h in ASW + Amicon + Mag_Beads", "Input_DNA", "Precipitation", "PCI + Precipitation", "Amicon + Precipitation", "PCI + Amicon + Precipitation"))
df.plot$Extraction <- factor(df.plot$Extraction, levels=c("Magnetic_Beads",
                                                              "Precipitation"))

df.plot$Extraction  <-  revalue(df.plot$Extraction ,c("Magnetic_Beads"="Magnetic_Beads"))
df$Extraction  <-  revalue(df$Extraction ,c("Magnetic_Beads"="Magnetic_Beads"))

# Calculate the "Percent" values based on "Value" and 2000 (maximum value)
df.plot$Percent <- df.plot$Value / 2000 * 100

# Plotting the data with secondary y-axis for "Percent" values

ggplot(df.plot, aes(Buffer, Value)) +
  geom_bar(stat = "identity", col = "black", position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymax = Value + SD, ymin = Value - SD), position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(data = df, aes(factor(Buffer), Value), position = position_dodge(width = 0.8), size = 1, color = "black") +
  geom_line(aes(y = Percent), linetype = "dashed", size = 1, color = "red") +
  facet_wrap(~ Extraction, scales = "free", ncol = 2) +
  labs(title = "Amount of recovered DNA - no mineral binding",
       y = "Total recovered DNA (ng)",
       ysec = "Percent of 2000") +  # Custom label for the secondary y-axis
  scale_y_continuous(
    limits = c(0, 2000),
    expand = c(0, 0.05),
    name = "Total recovered DNA (ng)",
    sec.axis = sec_axis(~ . / 20, name = "% input DNA recovered")) +
  theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.4, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(face = "bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size = 14),
        strip.background = element_blank(), strip.text = element_text(size = 12, margin = margin(2, 0, 2, 0)))
