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
library(ggpubr)

# Import the data
df <- read.table("smear_percentage.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

## melt the table into wide one
#wide_table <- dcast(df, Sample_ID ~ Range, value.var = "Value")
#write.table(wide, file = "wide_percentage_table", quote = FALSE, sep ="\t", row.names = FALSE)

# Drop failed samples
dropped <- c('405', '419', '420', '421')
df <- subset(df, !Sample_ID %in% dropped)

##### another way of testing is lme (more power) - we take sediment as a random effect

## Mixed linear model here gives us the option to use sediment_ID as a random effect as we do not have
## enough observations in a single sediment cluster. We thus use it as a random effect to increase observations for
## statistical comparison between the buffers while we disregard the method!

## Plot

ggplot(df, aes(x=Buffer, y=Value, group = Sediment_ID, color = Sediment_ID)) + 
  geom_smooth(method = lm, se = FALSE) +
  geom_point()+
  facet_grid( ~ Range) +
  theme_bw()

### Now we subset per range and test each range

df.1 <- df[df$Range == "157-177 bp",]
df.2 <- df[df$Range == "177-202 bp",]
df.3 <- df[df$Range == "202-227 bp",]
df.4 <- df[df$Range == "227-277 bp",]
df.5 <- df[df$Range == "277-327 bp",]
df.6 <- df[df$Range == "327-427 bp",]

### Build a model 1

model.1 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.1)
summary(model.1)
emmeans_model.1 <- emmeans(model.1, ~ Buffer)
print(emmeans_model.1)
posthoc.1 <- pairs(emmeans_model.1, adjust = "tukey")
print(posthoc.1)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   -0.147 0.0317 10.1  -4.634  0.0009

### Build a model 2

model.2 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.2)
summary(model.2)
emmeans_model.2 <- emmeans(model.2, ~ Buffer)
print(emmeans_model.2)
posthoc.2 <- pairs(emmeans_model.2, adjust = "tukey")
print(posthoc.2)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland  -0.0375 0.0352 10.4  -1.066  0.3105

### Build a model 3

model.3 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.3)
summary(model.3)
emmeans_model.3 <- emmeans(model.3, ~ Buffer)
print(emmeans_model.3)
posthoc.3 <- pairs(emmeans_model.3, adjust = "tukey")
print(posthoc.3)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0337 0.0182 10.2   1.856  0.0926

### Build a model 4

model.4 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.4)
summary(model.4)
emmeans_model.4 <- emmeans(model.4, ~ Buffer)
print(emmeans_model.4)
posthoc.4 <- pairs(emmeans_model.4, adjust = "tukey")
print(posthoc.4)

#contrast          estimate    SE   df t.ratio p.value
#Direito - Rohland   0.0615 0.018 10.1   3.421  0.0064

### Build a model 5

model.5 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.5)
summary(model.5)
emmeans_model.5 <- emmeans(model.5, ~ Buffer)
print(emmeans_model.5)
posthoc.5 <- pairs(emmeans_model.5, adjust = "tukey")
print(posthoc.5)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0418 0.0116 10.1   3.605  0.0047

### Build a model 6

model.6 <- lmer(Value ~ Buffer + (1 | Sediment_ID),
                data = df.6)
summary(model.6)
emmeans_model.6 <- emmeans(model.6, ~ Buffer)
print(emmeans_model.6)
posthoc.6 <- pairs(emmeans_model.6, adjust = "tukey")
print(posthoc.6)

#contrast          estimate     SE   df t.ratio p.value
#Direito - Rohland   0.0484 0.0162 10.1   2.995  0.0133


### four out of six ranges show significant difference between the two buffers in terms of molarity proportion of 
### different ranges

## Plot again and marking ranges as significant with a star

## Plot

plot <- ggplot(df, aes(x=Buffer, y=Value, group = Sediment_ID, color = Sediment_ID)) + 
  geom_smooth(method = lm, se = FALSE) +
  geom_point() +
  facet_grid( ~ Range) +
  theme_bw() + ylim(0,0.7) + ylab("Relative molarity [%]") + ggtitle("Relative molarity across fragment ranges")

plot

# Custom-defined text for each panel
panel_text <- data.frame(Range = c("157-177 bp", "177-202 bp", "202-227 bp", "227-277 bp", "277-327 bp", "327-427 bp"),
    label = c("p-value < 0.1", "ns", "ns", "p-value < 0.1", "p-value < 0.1", "p-value = 0.0133"))


plot + geom_text(data = panel_text, aes(x = 1.5, y = 0.7, label = label), vjust = 1.5, hjust = 0.5, inherit.aes = FALSE)


