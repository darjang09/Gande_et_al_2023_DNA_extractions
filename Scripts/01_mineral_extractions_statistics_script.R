## This is the script for statistical analysis of mineral-DNA extraction data and plotting of Figure 1

## Load packages

library(data.table)
library(ggplot2)
library(outliers)
library(car)
library(rcompanion)
library(FSA)
library(dplyr)

# Import the data

df <- read.table("data_for_stats_minerals_FINAL.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

# Inspect the data

ggplot(df, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + facet_grid(~Mineral)

boxplot(Value ~ Buffer, data = df, las = 2, xlab = "")

hist(df$Value)

############################### STATISTICS #####################################
################################################################################


########################## Scheirer–Ray–Hare Test ##############################
########################## non-parametric analysis #############################

## General combined data plot

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

# DV:  Value 
# Observations:  44 
# D:  0.9997181 
# MS total:  165 
# 
#                    Df Sum Sq     H p.value
# Buffer             4 5536.4 33.563 0.00000
# Extraction         1    9.7  0.059 0.80849
# Buffer:Extraction  4  297.1  1.801 0.77228
# Residuals         34 1255.0     

ggplot(ben, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

### KAOLIN ###

scheirerRayHare(Value ~ Buffer*Extraction, data = kao, type = 2)

# DV:  Value 
# Observations:  50 
# D:  0.9920768 
# MS total:  212.5 
# 
#                   Df Sum Sq      H p.value
# Buffer             4 8591.1 40.752 0.00000
# Extraction         1   21.8  0.103 0.74789
# Buffer:Extraction  4  276.2  1.310 0.85964
# Residuals         40 1440.9  

ggplot(kao, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

### SAND ###

scheirerRayHare(Value ~ Buffer*Extraction, data = san, type = 2)

# DV:  Value 
# Observations:  51 
# D:  1 
# MS total:  221 
# 
#                   Df Sum Sq       H p.value
# Buffer             4 5540.7 25.0711 0.00005
# Extraction         1   72.0  0.3260 0.56803
# Buffer:Extraction  4  247.5  1.1198 0.89112
# Residuals         41 5201.8   

ggplot(san, aes(Buffer, Value, fill=Extraction)) + geom_boxplot() + theme_bw()

### we see that there a significant effect of buffer for all three minerals and no effect of extraction type
### we can therefore follow up with a post-hoc test for the buffer

############################ POST-HOC TESTS ####################################
################################################################################

#____DUNN-TEST-BENTONITE_______________________________________________________#

## Order groups by median (here we use letters for buffers (see df) as R otherwise orders them alphabetically)

ben$Buffer = factor(ben$Buffer,
                    levels=c("A", "B", "C", "D","E"))

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

# Group Letter MonoLetter
# 1     A      a        a  
# 2     B      a        a  
# 3     C     ab        ab 
# 4     D     bc         bc
# 5     E      c          c

#____DUNN-TEST-KAOLIN__________________________________________________________#

## Order groups by median (here we use letters for buffers (see df) as R otherwise orders them alphabetically)

kao$Buffer = factor(kao$Buffer,
                    levels=c("A", "B", "C", "D","E"))

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

# Group Letter MonoLetter
# 1     A     ab        ab 
# 2     B      a        a  
# 3     C      c          c
# 4     D     bc         bc
# 5     E      c          c

#____DUNN-TEST-SAND____________________________________________________________#

## Order groups by median (here we use letters for buffers as R (see df) otherwise orders them alphabetically)

san$Buffer = factor(san$Buffer,
                    levels=c("A", "B", "C", "D","E"))

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

# Group Letter MonoLetter
# 1     A     ab         ab
# 2     B     abc         abc
# 3     C      c         c 
# 4     D      a          a
# 5     E      bc          bc

################################################################################
################################## PLOT CLD ####################################

#### BENTONITE

ben$Name <- factor(ben$Name, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display - here we now use buffer names instead of letters!

test_ben <- data.frame(Name=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("a", "a", "ab", "bc", "c"))

p1 <- ggplot(ben, aes(Name, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 7, aes(color=Extraction), show.legend = TRUE) +
  geom_text(data = test_ben, aes(label = labels, y = -110), show.legend = FALSE, size=14) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-140, 1000)) + 
  labs(title ="Bentonite", y = "Total recovered DNA [ng]") +
  theme(axis.text.x = element_text(margin = margin(t = 30, r = 0, b = 0.3, l = 0), size=34, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=34),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=34),
        legend.position="bottom",
        title =element_text(size=34, face='bold'),
        legend.text = element_text(size=34),
        legend.title = element_text(size=36),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 2))

p1 <- p1 + xlab("Buffer") + guides(color = guide_legend(title = "Isolation approach"))

#### KAOLIN

kao$Name <- factor(kao$Name, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display! - here we now use buffer names instead of letters!

test_kao <- data.frame(Name=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("ab", "a", "c", "bc", "c"))

p2 <- ggplot(kao, aes(Name, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 7, aes(color=Extraction), show.legend = TRUE) +
  geom_text(data = test_kao, aes(label = labels, y = -110), show.legend = FALSE, size=14) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-140, 1000)) + 
  labs(title ="Kaolinite", y = "") +
  theme(axis.text.x = element_text(margin = margin(t = 30, r = 0, b = 0.3, l = 0), size=34, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=36),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=36),
        legend.position="bottom",
        title =element_text(size=34, face='bold'),
        legend.text = element_text(size=34),
        legend.title = element_text(size=36),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 2))

p2 <- p2 + xlab("Buffer") + guides(color = guide_legend(title = "Isolation approach"))

#### SAND

san$Name <- factor(san$Name, levels=c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever"))

# assign labels according to compact letter display! - here we now use buffer names instead of letters!

test_san <- data.frame(Name=as.factor(c("Rohland", "Direito", "PowerSoil", "Pedersen", "Lever")), 
                       labels=c("ab", "abc", "c", "a", "bc"))

p3 <- ggplot(san, aes(Name, Value)) +
  geom_boxplot(position=position_dodge(0.8), outlier.alpha = 0, show.legend = FALSE) +
  geom_jitter(position=position_dodge(0.2), alpha = 0.5, cex = 7, aes(color=Extraction), show.legend = TRUE) +
  geom_text(data = test_san, aes(label = labels, y = -110), show.legend = FALSE, size=14) +
  theme_bw() +
  scale_color_manual(values=c("deepskyblue3", "orangered2")) +
  scale_y_continuous(expand=c(0,0), limits=c(-140, 1000)) + 
  labs(title ="Sea sand", y = "") +
  theme(axis.text.x = element_text(margin = margin(t = 30, r = 0, b = 0.3, l = 0), size=34, angle = 90, vjust=0.4, hjust=0, color = "black"),
        axis.text.y = element_text(size = 34, color = "black"),
        axis.title = element_text(size = 36, color = "black"),
        axis.title.x = element_text(margin = margin(t = 30, r = 0, b = 0, l = 0), size=36),
        axis.title.y = element_text(margin = margin(t = 0, r = 30, b = 0, l = 0),size=36),
        legend.position="bottom",
        title =element_text(size=34, face='bold'),
        legend.text = element_text(size=34),
        legend.title = element_text(size=36),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 2))

p3 <- p3 + xlab("Buffer") + guides(color = guide_legend(title = "Isolation approach"))

## JOIN ALL into a single graph

library(ggpubr)

figure <- ggarrange(p1, p3, p2,
                    labels = c("A", "B", "C"),
                    font.label = list(size = 36, color = "black"),
                    ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom", align = "v")
figure
