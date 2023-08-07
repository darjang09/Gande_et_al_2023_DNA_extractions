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

qPCR <- read.table("../../All_the_data/06_Results/01_extraction_optimization/Extractions_2022/Statistics/data_for_stats_controls_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

agg=aggregate(Value~Buffer*Extraction, data=qPCR, FUN="mean") #mean
agg$sd=aggregate(Value~Buffer*Extraction, data=qPCR, FUN="sd")$Value #add the SD
agg$se <- with(agg, sd/2) #add the standard error (SD/2) #2 bc we have a root of 4 samples

qPCR.plot <- agg

names(qPCR.plot)[c(1:4)] <- c("Buffer","Extraction","Value","SD")
str(qPCR.plot)
qPCR.plot$Buffer <- factor(qPCR.plot$Buffer, levels=c("Input_DNA", "Mag_Beads_only", "Amicon + Mag_Beads", "24h in ASW + Amicon + Mag_Beads",
                                                      "Input_DNA", "Precipitation", "PCI + Precipitation", "Amicon + Precipitation", "PCI + Amicon + Precipitation"),
                           labels=c("Input_DNA", "Mag_Beads_only", "Amicon + Mag_Beads", "24h in ASW + Amicon + Mag_Beads", "Input_DNA", "Precipitation", "PCI + Precipitation", "Amicon + Precipitation", "PCI + Amicon + Precipitation"))
qPCR.plot$Extraction <- factor(qPCR.plot$Extraction, levels=c("Magnetic_Beads",
                                                              "Precipitation"))

qPCR.plot$Extraction  <-  revalue(qPCR.plot$Extraction ,c("Magnetic_Beads"="Magnetic_Beads"))
qPCR$Extraction  <-  revalue(qPCR$Extraction ,c("Magnetic_Beads"="Magnetic_Beads"))

# ggplot(qPCR.plot, aes(Buffer, Value)) +
#   geom_bar(stat="identity", col="black", position = position_dodge(width=0.8), width=0.8) +
#   scale_y_continuous(expand=c(0,0), limits=c(0, 2000)) +
#   geom_errorbar(aes(ymax = Value + SD, ymin = Value - SD), position = position_dodge(width=0.8), width=0.2) +
#   geom_point(data=qPCR, aes(factor(Buffer), Value), position = position_dodge(width=0.8), size=1, color="black") +
#   scale_fill_manual("Extraction method") +
#   facet_wrap( ~  Extraction , scales = "free", ncol=4) +
#   labs(title ="Amount of recovered DNA - no mineral binding",
#        y = "Total recovered DNA (ng)") +
#   theme(axis.text.x = element_text(size=11, angle = 90, vjust=0.4, hjust=1, color = "black"),
#         axis.text.y = element_text(size = 11, color = "black"),
#         axis.title = element_text(face="bold", size = 12, color = "black"),
#         axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
#         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(colour = "black", fill = "NA"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.text = element_text(size=12),
#         legend.title = element_text(size=14),
#         legend.key = element_rect(colour = NA, fill = NA),
#         legend.key.size = unit(0.5, "cm"),
#         plot.title = element_text(size=14),
#         strip.background = element_blank(), strip.text = element_text( size = 12, margin = margin(2,0,2,0)))
# 

##### ChatGPT created code to update the graph with secondary axis to plot percentage values

# Calculate the "Percent" values based on "Value" and 2000 (maximum value)
qPCR.plot$Percent <- qPCR.plot$Value / 2000 * 100

# Plotting the data with secondary y-axis for "Percent" values
ggplot(qPCR.plot, aes(Buffer, Value)) +
  geom_bar(stat = "identity", col = "black", position = position_dodge(width = 0.8), width = 0.8) +
  geom_errorbar(aes(ymax = Value + SD, ymin = Value - SD), position = position_dodge(width = 0.8), width = 0.2) +
  geom_point(data = qPCR, aes(factor(Buffer), Value), position = position_dodge(width = 0.8), size = 1, color = "black") +
  geom_line(aes(y = Percent), linetype = "dashed", size = 1, color = "red") +
  facet_wrap(~ Extraction, scales = "free", ncol = 2) +
  labs(title = "Amount of recovered DNA - no mineral binding",
       y = "Total recovered DNA (ng)",
       ysec = "Percent of 2000") +  # Custom label for the secondary y-axis
  scale_y_continuous(
    limits = c(0, 2000),
    expand = c(0, 0.05),
    name = "Total recovered DNA (ng)",
    sec.axis = sec_axis(~ . / 20, name = "% input DNA recovered")
  ) +
  
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

#table <- round(agg$sd, 2)

