## This is the final script for statistical analysis of mineral-DNA extraction data

## Load packages

library(data.table)
library(ggplot2)
library(outliers)
library(car)
library(rcompanion)
library(FSA)
library(dplyr)

# Import the data

df <- read.table("D:/All_the_data/All_the_data/06_Results/01_extraction_optimization/Extractions_2022/Statistics/data_for_stats_minerals_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

# Inspect the data

ggplot(df, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + facet_grid(~Mineral)

boxplot(Value ~ Buffer, data = df, las = 2, xlab = "")

hist(df$Value)

# We see that our data has non-normal distribution (even when values are transformed to log) which makes it hard for ANOVA
# to be reliable (not shown in this script). In addition our design is unbalanced!

## We therefore continue with other statistical tools to test our hypotheses
#______________________________________________________________________________#

########################### OUTLIER DETECTION ##################################

## Here we use Dixon test to identify potential outliers that would have a big impact on our data

# Overview
df$group <- as.factor(paste(df$Buffer, df$Extraction, df$Mineral, sep="_"))
table(df$group)

### TEST ###

dixon.test(df[df$group == "PowerSoil_Precipitation_Bentonite",]$Value, opposite = F)  #p < 0.05, most extreme value is an outlier

#the test takes the most extreme value into account and not the highest.
#With option opposite=F you specifiy the opposite of what it chooses by default

# If there is an outlier, remove it and test again
#v <- df[df$group == "PowerSoil_Precipitation_Bentonite",]$Value
#v <- v[v != max(v)]                                                            #in both cases our outlier was highest value
#dixon.test(v)

### RESULTS ###

#Rohland_Magnetic_Beads_Bentonite - p-value = 0.7064
#Rohland_Magnetic_Beads_Sand - p-value = 0.5237
#Rohland_Magnetic_Beads_Kaolin - p-value = 0.9795
#Rohland_Precipitation_Bentonite - p-value = 0.9208
#Rohland_Precipitation_Sand - p-value = 0.6552
#Rohland_Precipitation_Kaolin - p-value = 0.4766
#PowerSoil_Magnetic_Beads_Bentonite - p-value = 0.01189 ........................<- -max = p-value = 0.1384
#PowerSoil_Magnetic_Beads_Sand - p-value = 0.2025
#PowerSoil_Magnetic_Beads_Kaolin - cannot be calculated = 0
#PowerSoil_Precipitation_Bentonite - p-value = 0.01615 .........................<- -max = p-value = 0.1686
#PowerSoil_Precipitation_Sand - p-value = 0.3158
#PowerSoil_Precipitation_Kaolin - p-value = 0.7836
#Pedersen_Magnetic_Beads_Bentonite - p-value = 0.279
#Pedersen_Magnetic_Beads_Sand - p-value = 0.556
#Pedersen_Magnetic_Beads_Kaolin - p-value = 0.9179
#Pedersen_Precipitation_Bentonite - p-value = 0.1548
#Pedersen_Precipitation_Sand - p-value = 0.2841
#Pedersen_Precipitation_Kaolin - p-value = 0.3158
#Lever_Magnetic_Beads_Bentonite - cannot be calculated = 0
#Lever_Magnetic_Beads_Sand - p-value = 0.3184
#Lever_Magnetic_Beads_Kaolin - p-value = 0.4174
#Lever_Precipitation_Bentonite - p-value = 0.5299
#Lever_Precipitation_Sand - p-value = 0.1336
#Lever_Precipitation_Kaolin - p-value = 0.7395
#Direito_Magnetic_Beads_Bentonite -  p-value = 0.9825
#Direito_Magnetic_Beads_Sand - p-value = 0.1129
#Direito_Magnetic_Beads_Kaolin - p-value = 0.9172
#Direito_Precipitation_Bentonite - p-value = 0.8398
#Direito_Precipitation_Sand - p-value = 0.6432
#Direito_Precipitation_Kaolin - p-value = 0.1537

#______________________________________________________________________________#
# We can choose to remove or keep the two outliers in the data table before proceeding

df <- df[-c(27,42),]

############################### STATISTICS #####################################
################################################################################


# outliers removed (y/n) -> YES

########################## Scheirer–Ray–Hare Test ##############################
########################## non-parametric analysis #############################

## General combined data plot

df$Buffer <- factor(df$Buffer, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"),
                    labels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

df$Mineral <- factor(df$Mineral, levels=c("Bentonite", "Kaolin", "Sand"))

df$Extraction <- factor(df$Extraction, levels=c("Magnetic_Beads", "Precipitation"))

ggplot(df, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + facet_grid(~Mineral) +
  theme_bw() + scale_y_continuous(expand=c(0,0), limits=c(0, 2000)) + 
  labs(title ="Amount of recovered DNA", y = "Total recovered DNA [ng]") +
  theme(axis.text.x = element_text(size=11, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(face="bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12),
        panel.background = element_rect(fill = "white"))


## Sheier-Ray-Hare is only a two way test so we break the table per mineral

ben <- df[df$Mineral == 'Bentonite', ]
san <- df[df$Mineral == 'Sand', ]
kao <- df[df$Mineral == 'Kaolin', ]

### BENTONITE ###

scheirerRayHare(Value ~ Buffer*Extraction, data = ben, type = 2)

# WITH OUTLIERS
# DV:  Value 
# Observations:  44 
# D:  0.9985906 
# MS total:  165 
# 
#                   Df Sum Sq      H p.value
# Buffer             4 5390.0 32.713 0.00000
# Extraction         1   32.1  0.195 0.65900
# Buffer:Extraction  4  258.9  1.571 0.81392
# Residuals         34 1409.6 

# WITHOUT OUTLIERS
# DV:  Value 
# Observations:  42 
# D:  0.9983794 
# MS total:  150.5 
# 
# Df Sum Sq      H p.value
# Buffer             4 4981.4 33.153 0.00000
# Extraction         1   62.4  0.415 0.51937
# Buffer:Extraction  4  177.8  1.184 0.88079
# Residuals         32  933.8

ggplot(ben, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

### KAOLIN ###

scheirerRayHare(Value ~ Buffer*Extraction, data = kao, type = 2)

# DV:  Value 
# Observations:  50 
# D:  0.9983193 
# MS total:  212.5 
# 
#                   Df Sum Sq      H p.value
# Buffer             4 8954.1 42.208 0.00000
# Extraction         1   79.4  0.374 0.54073
# Buffer:Extraction  4  232.0  1.094 0.89528
# Residuals         40 1129.6    

ggplot(kao, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

### SAND ###

scheirerRayHare(Value ~ Buffer*Extraction, data = san, type = 2)

# DV:  Value 
# Observations:  51 
# D:  1 
# MS total:  221 
# 
#                    Df Sum Sq      H p.value
# Buffer             4 4300.2 19.4578 0.00064
# Extraction         1  218.1  0.9869 0.32049
# Buffer:Extraction  4 2042.5  9.2420 0.05533
# Residuals         41 4180.9

ggplot(san, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

#______________________________________________________________________________#

### we see that there a significant effect of buffer for all three minerals and no effect of extraction type
### we can therefore follow up with a post-hoc test for the buffer

############################ POST-HOC TESTS ####################################
################################################################################


#____DUNN-TEST-BENTONITE_______________________________________________________#

## Order groups by median

ben$Buffer = factor(ben$Buffer,
                    levels=c("Rohland", "Pedersen", "PowerSoil", "Lever","Direito"))

levels(ben$Buffer)

## Test

dt = dunnTest(Value ~ Buffer,
              data=ben,
              method="holm")      # here we choose "holm" to adjust p-values for multiple comparisons
dt

## Compact letter display

pt = dt$res

cldList(P.adj ~ Comparison,
        data = pt,
        threshold = 0.05)

# WITH OUTLIERS
#     Group   Letter MonoLetter
# 1   Direito      a        a  
# 2     Lever      b         b 
# 3  Pedersen     bc         bc
# 4 PowerSoil     ac        a c
# 5   Rohland      a        a 

# WITHOUT OUTLIERS
#       Group Letter MonoLetter
# 1   Direito      a         a 
# 2     Lever      b          b
# 3  Pedersen      b          b
# 4 PowerSoil     ab         ab
# 5   Rohland      a         a

#____DUNN-TEST-KAOLIN__________________________________________________________#

## Order groups by median

kao$Buffer = factor(kao$Buffer,
                    levels=c("Rohland", "Pedersen", "PowerSoil", "Lever","Direito"))

levels(kao$Buffer)

## Test

dt = dunnTest(Value ~ Buffer,
              data=kao,
              method="holm")      # here we choose "holm" to adjust p-values for multiple comparisons
dt

## Compact letter display

pt = dt$res

cldList(P.adj ~ Comparison,
        data = pt,
        threshold = 0.05)

#       Group Letter MonoLetter
# 1   Direito      a        a  
# 2     Lever      b         b 
# 3  Pedersen     bc         bc
# 4 PowerSoil      b         b 
# 5   Rohland     ac        a c

#____DUNN-TEST-SAND____________________________________________________________#

san$Buffer = factor(san$Buffer,
                    levels=c("Rohland", "Pedersen", "PowerSoil", "Lever","Direito"))

levels(san$Buffer)

## Test

dt = dunnTest(Value ~ Buffer,
              data=san,
              method="holm")      # here we choose "holm" to adjust p-values for multiple comparisons
dt

## Compact letter display

pt = dt$res

cldList(P.adj ~ Comparison,
        data = pt,
        threshold = 0.05)

#       Group Letter MonoLetter
# 1   Direito     ab         ab
# 2     Lever      a         a 
# 3  Pedersen      a         a 
# 4 PowerSoil      b          b
# 5   Rohland     ab         ab

################################################################################
################################## PLOT CLD ####################################

#### BENTONITE

# custom order buffers

ben$Buffer <- factor(ben$Buffer, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display!

test_ben <- data.frame(Buffer=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("a", "a", "ab", "b", "b")) # you can extract the table from actual dunn test!

p1 <- ggplot(ben, aes(Buffer, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 3, aes(color=Extraction), show.legend = TRUE) +
  geom_text(data = test_ben, aes(label = labels, y = -75), show.legend = FALSE) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-150, 1500)) + 
  labs(title ="Bentonite", y = "Total recovered DNA [ng]") +
  theme(axis.text.x = element_text(size=11, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(face="bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12),
        legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"))
p1
#### KAOLIN

# custom order buffers

kao$Buffer <- factor(kao$Buffer, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display!

test_kao <- data.frame(Buffer=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("ac", "a", "b", "bc", "b")) # you can extract the table from actual dunn test!

p2 <- ggplot(kao, aes(Buffer, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 3, aes(color=Extraction), show.legend = TRUE) +
  geom_text(data = test_kao, aes(label = labels, y = -75), show.legend = FALSE) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-150, 2000)) + 
  labs(title ="Kaolinite", y = "") +
  theme(axis.text.x = element_text(size=11, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(face="bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12),
        legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"))
p2
#### SAND

# custom order buffers

san$Buffer <- factor(san$Buffer, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display!

test_san <- data.frame(Buffer=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("ab", "ab", "b", "a", "a")) # you can extract the table from actual dunn test!

p3 <- ggplot(san, aes(Buffer, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 3, aes(colour=Extraction),  show.legend = TRUE) +
  geom_text(data = test_san, aes(label = labels, y = -75), show.legend = FALSE) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-150, 2000)) + 
  labs(title ="Sand", y = "") +
  theme(axis.text.x = element_text(size=11, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title = element_text(face="bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=12),
        legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"))
p3

## JOIN ALL into a single graph

library(ggpubr)

figure <- ggarrange(p1, p3, p2,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom")
figure