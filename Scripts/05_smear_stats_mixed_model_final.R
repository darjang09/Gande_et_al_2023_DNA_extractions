## This script was used for in depth statistical analysis of fragment analyzer data and additional plotting of Figure 4 and Figure S6

## Load libraries

library("data.table")
library("ggplot2")
library("car")
library("FSA")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("multcomp")
library("sandwich")
library("reshape2")
library("utils")
library("lme4")
library("emmeans")
library("ggpubr")

# Import the data
df <- read.table("smear_percentage_and_absolute.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

# Drop failed samples
dropped <- c('405', '419', '420', '421')
df <- subset(df, !Sample_ID %in% dropped)

################################################################################
############# Linear Mixed Effects Models With Relative Values #################

### We first subset the data per range and test each range separately

df.1 <- df[df$Range == "157-177 bp",]
df.2 <- df[df$Range == "177-202 bp",]
df.3 <- df[df$Range == "202-227 bp",]
df.4 <- df[df$Range == "227-277 bp",]
df.5 <- df[df$Range == "277-327 bp",]
df.6 <- df[df$Range == "327-427 bp",]

### Build model 1

model.1 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.1)
summary(model.1)
emmeans_model.1 <- emmeans(model.1, ~ Buffer)
print(emmeans_model.1)
posthoc.1 <- pairs(emmeans_model.1, adjust = "tukey")
print(posthoc.1)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   -0.147 0.0317 10.1  -4.634  0.0009 <------<------<------

### Build model 2

model.2 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.2)
summary(model.2)
emmeans_model.2 <- emmeans(model.2, ~ Buffer)
print(emmeans_model.2)
posthoc.2 <- pairs(emmeans_model.2, adjust = "tukey")
print(posthoc.2)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland  -0.0375 0.0352 10.4  -1.066  0.3105

### Build model 3

model.3 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.3)
summary(model.3)
emmeans_model.3 <- emmeans(model.3, ~ Buffer)
print(emmeans_model.3)
posthoc.3 <- pairs(emmeans_model.3, adjust = "tukey")
print(posthoc.3)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0337 0.0182 10.2   1.856  0.0926

### Build model 4

model.4 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.4)
summary(model.4)
emmeans_model.4 <- emmeans(model.4, ~ Buffer)
print(emmeans_model.4)
posthoc.4 <- pairs(emmeans_model.4, adjust = "tukey")
print(posthoc.4)

#contrast          estimate    SE   df t.ratio p.value
#Direito - Rohland   0.0615 0.018 10.1   3.421  0.0064 <------<------<-------

### Build model 5

model.5 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.5)
summary(model.5)
emmeans_model.5 <- emmeans(model.5, ~ Buffer)
print(emmeans_model.5)
posthoc.5 <- pairs(emmeans_model.5, adjust = "tukey")
print(posthoc.5)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0418 0.0116 10.1   3.605  0.0047 <------<------<------

### Build a model 6

model.6 <- lmer(Percent ~ Buffer + (1 | Sediment_ID),
                data = df.6)
summary(model.6)
emmeans_model.6 <- emmeans(model.6, ~ Buffer)
print(emmeans_model.6)
posthoc.6 <- pairs(emmeans_model.6, adjust = "tukey")
print(posthoc.6)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0484 0.0162 10.1   2.995  0.0133 <-------<-------<-------

### Four out of six ranges show significant difference between the two buffers in terms of relative molarity

### Lets plot and mark significant ranges

p1 <- ggplot(df, aes(x=Buffer, y=Percent, group = Sediment_ID, color = Sediment_ID)) + 
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  facet_grid( ~ Range) +
  theme_minimal() +
  ylim(0,0.7) + ylab("Relative molarity [%]") +
  ggtitle("Relative molarity across fragment size intervals") +
  theme(legend.position="bottom",
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 0.3, l = 0), size=12, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=14),
        title =element_text(size=14, face='bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

p1

# Add p-values for each panel
panel_text <- data.frame(Range = c("157-177 bp", "177-202 bp", "202-227 bp", "227-277 bp", "277-327 bp", "327-427 bp"),
    label = c("p < 0.01", "ns", "ns", "p < 0.01", "p < 0.01", "p = 0.0133"))

p1 <- p1 + geom_text(data = panel_text, aes(x = 1.5, y = 0.7, label = label), vjust = 1.5, hjust = 0.5, inherit.aes = FALSE)

p1
################################################################################
############# Linear Mixed Effects Models With Absolute Values #################

### Build model 1

model.1 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.1)
summary(model.1)
emmeans_model.1 <- emmeans(model.1, ~ Buffer)
print(emmeans_model.1)
posthoc.1 <- pairs(emmeans_model.1, adjust = "tukey")
print(posthoc.1)

#contrast          estimate   SE   df t.ratio p.value
#Direito - Rohland   -0.409 2.14 10.3  -0.191  0.8525

### Build model 2

model.2 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.2)
summary(model.2)
emmeans_model.2 <- emmeans(model.2, ~ Buffer)
print(emmeans_model.2)
posthoc.2 <- pairs(emmeans_model.2, adjust = "tukey")
print(posthoc.2)

#contrast          estimate   SE   df t.ratio p.value
#Direito - Rohland      4.3 3.97 11.4   1.084  0.3009

### Build model 3

model.3 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.3)
summary(model.3)
emmeans_model.3 <- emmeans(model.3, ~ Buffer)
print(emmeans_model.3)
posthoc.3 <- pairs(emmeans_model.3, adjust = "tukey")
print(posthoc.3)

#contrast          estimate   SE   df t.ratio p.value
#Direito - Rohland     2.77 1.85 11.3   1.493  0.1628

### Build a model 4

model.4 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.4)
summary(model.4)
emmeans_model.4 <- emmeans(model.4, ~ Buffer)
print(emmeans_model.4)
posthoc.4 <- pairs(emmeans_model.4, adjust = "tukey")
print(posthoc.4)

#contrast          estimate   SE   df t.ratio p.value
#Direito - Rohland     2.83 1.67 10.8   1.694  0.1189

### Build a model 5

model.5 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.5)
summary(model.5)
emmeans_model.5 <- emmeans(model.5, ~ Buffer)
print(emmeans_model.5)
posthoc.5 <- pairs(emmeans_model.5, adjust = "tukey")
print(posthoc.5)

#contrast          estimate   SE   df t.ratio p.value
#Direito - Rohland     1.46 0.68 10.6   2.150  0.0555

### Build a model 6

model.6 <- lmer(Absolute ~ Buffer + (1 | Sediment_ID),
                data = df.6)
summary(model.6)
emmeans_model.6 <- emmeans(model.6, ~ Buffer)
print(emmeans_model.6)
posthoc.6 <- pairs(emmeans_model.6, adjust = "tukey")
print(posthoc.6)

#contrast          estimate    SE   df t.ratio p.value
#Direito - Rohland     1.46 0.564 10.4   2.584  0.0264 <-------<-------<--------   

### one out of six ranges show significant difference between the two buffers in terms of absoulte molarity different

## Plot

p2 <- ggplot(df, aes(x=Buffer, y=Absolute, group = Sediment_ID, color = Sediment_ID)) + 
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  facet_grid( ~ Range) +
  theme_minimal() +
  ylim(0,40) + ylab("Absolute molarity [nM]") +
  ggtitle("Absolute molarity across fragment size intervals") +
  theme(legend.position="bottom",
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 0.3, l = 0), size=12, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=14),
        title =element_text(size=14, face='bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

p2

# Add p-values for each panel
panel_text <- data.frame(Range = c("157-177 bp", "177-202 bp", "202-227 bp", "227-277 bp", "277-327 bp", "327-427 bp"),
                         label = c("ns", "ns", "ns", "ns", "ns", "p = 0.0264"))


p2 <- p2 + geom_text(data = panel_text, aes(x = 1.5, y = 40, label = label), vjust = 1.5, hjust = 0.5, inherit.aes = FALSE)

p2

### Merge both plots for Figure S6

figure <- ggarrange(p1, p2,
                    labels = c("A", "B"),
                    font.label = list(size = 20, color = "black", face = "bold"),
                    ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom", align = "v", hjust = 0, vjust = 1.2)
figure

################################################################################
################ Plotting Figure 4 for indepth interpretation ##################

p3 <- ggplot(df, aes(Method, Percent, color = Buffer)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.8), alpha = 0.5, cex = 4, aes(colour=Buffer), show.legend = TRUE) +
  facet_grid( ~ Range) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous() + 
  labs(title ="Relative molarity across fragment size intervals", y = "Relative molarity [%]", x = "Isolation approach") +
  theme(legend.position="bottom",
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 0.3, l = 0), size=12, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=14),
        title =element_text(size=14, face='bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))

p4 <- ggplot(df, aes(Method, Absolute, color = Buffer)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.8), alpha = 0.5, cex = 4, aes(colour=Buffer), show.legend = TRUE) +
  facet_grid( ~ Range) +
  theme_minimal() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous() + 
  labs(title ="Molarity across fragment size intervals", y = "Molarity [nM]", x = "Isolation approach") +
  theme(legend.position="bottom",
        axis.ticks.y = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(margin = margin(t = 0, r = 0, b = 0.3, l = 0), size=12, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=14),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=14),
        title =element_text(size=14, face='bold'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14))


figure1 <- ggarrange(p3, p4,
                 labels = c("A", "B"),
                 font.label = list(size = 20, color = "black", face = "bold"),
                 ncol = 1, nrow = 2, common.legend = TRUE, legend="bottom", align = "v", hjust = 0, vjust = 1.2)
figure1
