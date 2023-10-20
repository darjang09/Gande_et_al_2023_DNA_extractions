## Script for analysis and plotting of fragment analyzer data - Figure 3

## Load libraries

library(data.table)
library(XML)
library(base64enc)
library(ggplot2)
library(tidyverse)
library(imputeTS)
library(reshape2)
library(dplyr)
library(ggpubr)
library(scales)
library(graphics)

## Load the data
df <- read.table("merged_smear_April23.csv", sep=";", header=T, stringsAsFactors = T, check.names = FALSE)

## Rename range names
df <- df %>% mutate(Range = gsub(" bp to ", "-", Range))

############################## PARSE TABLE #####################################

### Subset columns
df <- df %>% select(`Sample ID`, Range, `nmole/L`) # for molarity

### Melt the table
data.perc <- reshape(df, idvar = "Range", timevar = "Sample ID", direction = "wide")
data.perc <- data.perc %>% remove_rownames %>% column_to_rownames(var="Range")

## Here we exclude 25-157 bp and 157-500 bp from the data since we don´t need it in this analysis
data.perc <- data.perc[row.names(data.perc) != "25-157 bp", , drop = FALSE]
data.perc <- data.perc[row.names(data.perc) != "157-500 bp", , drop = FALSE]

## Trim names
names(data.perc) = gsub(pattern = "*nmole/L.", replacement = "", x = names(data.perc)) # for molarity

## Transform values back to numerical
data.perc <- mutate_all(data.perc, function(x) as.numeric(as.character(x)))

## Transform NaN values to 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data.perc[is.nan(data.perc)] <- 0

## Make a table for absolute molarity
data.absolute <- as.matrix(data.perc)
data.absolute <- melt(data.absolute, id.vars=c("bracket"))

### FINALLY we transform everything into a percentage
data.perc <- apply(data.perc, 2, function(x){x/sum(x)})
colSums(data.perc) # check if all values are 1!
data.perc <- melt(data.perc, id.vars=c("bracket"))

##### Define sediment type list for faceting the barplots later

sedlist <- list(
  N_PAC = c('400', '412', '405','417'),
  MID_PAC = c('401','413','406','418'),
  BS_SAP = c('402','414','407','419'),
  BS_LAC = c('403','415','408','420'),
  SG_CLY = c('404','416','409','421'),
  EXT_BLANK = c('410','422','411','423'),
  LIB_BLANK = c('Lib_b_1','Lib_b_2','Lib_b_3','Lib_b_4')
)

##################### Plotting relative molarity #########################

## define color blind friendly pallete

cols <- palette(hcl.colors(6, "heat"))

level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','422','411','423','Lib_b_1','Lib_b_2','Lib_b_3','Lib_b_4')

sedlevel <- c('N_PAC', 'MID_PAC', 'BS_SAP', 'BS_LAC', 'SG_CLY', 'EXT_BLANK', 'LIB_BLANK')

data.perc$Sediment <- factor(data.perc$Var2)
levels(data.perc$Sediment) <- sedlist

p1 <- ggplot(data.perc, aes(fill=Var1, y=value, x=factor(Var2, level = level_order))) +
  geom_bar(position="fill", stat="identity", col="Black") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size=28, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title = element_text(face="bold", size = 30, color = "black"),
        #axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=28),
        legend.text = element_text(size=28),
        legend.title = element_text(size=30),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=28),
        strip.text.x = element_text(size = 28)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = FALSE),limits = force) +
  ylab("Relative Molarity [%]") + xlab("Sample ID") +
  guides(fill=guide_legend(title="Fragment range")) +
  facet_grid(. ~ factor(Sediment, level = sedlevel), scales = "free", space = "free")

p1

##################### Plotting absolute molarity ###############################

get_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n")
  }
}

level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','422','411','423','Lib_b_1','Lib_b_2','Lib_b_3','Lib_b_4')

cols <- palette(hcl.colors(6, "heat"))

sedlevel <- c('N_PAC', 'MID_PAC', 'BS_SAP', 'BS_LAC', 'SG_CLY', 'EXT_BLANK', 'LIB_BLANK')

data.absolute$Sediment <- factor(data.absolute$Var2)
levels(data.absolute$Sediment) <- sedlist
data.absolute$Var1 <- factor(data.absolute$Var1, levels=c('157-177 bp', '177-202 bp', '202-227 bp', '227-277 bp', '277-327 bp', '327-427 bp'))

p2 <- ggplot(data.absolute, aes(x=factor(Var2, level = level_order), y=value, fill=Var1)) + 
  geom_bar(position='stack', stat='identity', col="Black") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size=28, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title = element_text(face="bold", size = 30, color = "black"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size=28),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=28),
        legend.text = element_text(size=28),
        legend.title = element_text(size=30),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=28),
        strip.text.x = element_text(size = 28)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,65)) +
  scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = FALSE),limits = force) +
  ylab("Molarity [nM]") + xlab("") +
  guides(fill=guide_legend(title="Fragment range")) +
  facet_grid(. ~ factor(Sediment, level = sedlevel),  scales = "free", space = "free")

p2
                        
#########################################################################################
############################ PLOTTING THE RATIO #########################################

ratio <- read.table("FA_molarity_ratios_25-157-500.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

## facet them per sediment type

sedlevel <- c('N_PAC', 'MID_PAC', 'BS_SAP', 'BS_LAC', 'SG_CLY', 'EXT_BLANK', 'LIB_BLANK')
level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','422','411','423','Lib_b_1','Lib_b_2','Lib_b_3','Lib_b_4') 

p3 <- ggplot(ratio, aes(x=factor(name, level = level_order), y=ratio, shape=Extraction, fill=Buffer, group=name)) +
  geom_point(size = 6) + ylim(0,4) +
  geom_linerange(aes(x=factor(name, level = level_order), ymax=ratio, ymin=0, group=name), position = "identity") +
  scale_shape_manual(values = c(21, 22, 23), breaks=c('Magnetic Beads', 'Precipitation', 'Blank')) +
  scale_fill_manual(
    values = c("red", "blue", "black"),
    breaks=c('Rohland', 'Direito', 'Blank'),
    guide = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.x = element_text(size=28, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title = element_text(face="bold", size = 30, color = "black"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0), size=28),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0), size=28),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        legend.position = "right",
        legend.text = element_text(size=28),
        legend.title = element_text(size=30),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        strip.background = element_blank(), strip.text = element_text(size = 28, margin = margin(2,0,2,0))) +
  labs(x ="", y = "Molarity ratio\n[157-500 bp] / [25-157 bp]") +
  facet_grid(. ~ factor(Sediment, level = sedlevel), scales = "free", space = "free")

p3

#### MERGE ALL GRAPHS ####
library(ggpubr)

figure <- ggarrange(p3, p2, p1,
                    labels = c("A", "B", "C"),
                    font.label = list(size = 32, color = "black", face = "bold"),
                    ncol = 1, nrow = 3, common.legend = FALSE, legend="right", align = "v", hjust = 0, vjust = 1.2)
figure
