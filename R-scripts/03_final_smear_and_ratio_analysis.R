## Load libraries
library(bioanalyzeR)
library(data.table)
library(XML)
library(base64enc)
library(ggplot2)
library(tidyverse)
library(imputeTS)
library(reshape2)
library(dplyr)

## Load the data

df <- read.table("D:/All_the_data/Current_R_projects/EXTRACTIONS_ANALYSIS/merged_smear_April23.csv", sep=";", header=T, stringsAsFactors = T, check.names = FALSE)

############ PARSE TABLE

### subset columns , melt the table, tidy up namings and transform into %
df <- df %>% select(`Sample ID`, Range, `nmole/L`)
#df <- df %>% select(`Sample ID`, Range, `ng/uL`)
data.perc <- reshape(df, idvar = "Range", timevar = "Sample ID", direction = "wide")
data.perc <- data.perc %>% remove_rownames %>% column_to_rownames(var="Range")

## we exclude 25-157 bp and 157-500 bp from the data since it´s useless for us to determine molarity of sequencable fragments

data.perc <- data.perc[row.names(data.perc) != "25 bp to 157 bp", , drop = FALSE]
data.perc <- data.perc[row.names(data.perc) != "157 bp to 500 bp", , drop = FALSE]

## trim names
names(data.perc) = gsub(pattern = "*nmole/L.", replacement = "", x = names(data.perc))
#names(data.perc) = gsub(pattern = "*ng/uL.", replacement = "", x = names(data.perc))

## transform values back to numerical
data.perc <- mutate_all(data.perc, function(x) as.numeric(as.character(x)))

## transform NaN values to 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data.perc[is.nan(data.perc)] <- 0


## make a table for absolute molarity
data.absolute <- as.matrix(data.perc)
data.absolute <- melt(data.absolute, id.vars=c("bracket"))

### FINALLY transform everything into a percentage
data.perc <- apply(data.perc, 2, function(x){x/sum(x)})
colSums(data.perc)
data.perc <- melt(data.perc, id.vars=c("bracket"))


##### define sediment type list for faceting the barplots later

sedlist <- list(
  N_PAC = c('400', '412', '405','417'),
  MID_PAC = c('401','413','406','418'),
  BS_SAP = c('402','414','407','419'),
  BS_LAC = c('403','415','408','420'),
  SG_CLY = c('404','416','409','421'),
  EXTRACTION_BLANK = c('410','411','422','423'),
  LIBRARY_BLANK = c('Lib_blank_1(1-12)','Lib_blank_2(1-12)','Lib_blank_1(13-24)','Lib_blank_2(13-24)')
)

##################### Plotting percentage molarity brackets

cols <- c("#005db0", "#84005f", "#bd0040", "#eb5248", "#f48734", "#f0d35b")
level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','411','422','423','Lib_blank_1(1-12)','Lib_blank_2(1-12)','Lib_blank_1(13-24)','Lib_blank_2(13-24)')
sedlevel <- c('N_PAC', 'MID_PAC', 'BS_SAP', 'BS_LAC', 'SG_CLY', 'EXTRACTION_BLANK', 'LIBRARY_BLANK')

data.perc$Sediment <- factor(data.perc$Var2)
levels(data.perc$Sediment) <- sedlist

ggplot(data.perc, aes(fill=Var1, y=value, x=factor(Var2, level = level_order))) +
  geom_bar(position="fill", stat="identity", col="Black") +
  theme_bw() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(face="bold", size = 18, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=18),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=14),
        strip.text.x = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = FALSE),limits = force) +
  ylab("Molarity (%)") + xlab("Sample") +
  guides(fill=guide_legend(title="Fragment length")) +
  facet_grid(. ~ factor(Sediment, level = sedlevel), scales = "free", space = "free")


##################### Plotting absolute molarity

get_wraper <- function(width) {
  function(x) {
    lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n")
  }
}

level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','411','422','423','Lib_blank_1(1-12)','Lib_blank_2(1-12)','Lib_blank_1(13-24)','Lib_blank_2(13-24)')
cols <- c("#005db0", "#84005f", "#bd0040", "#eb5248", "#f48734", "#f0d35b")
sedlevel <- c('N_PAC', 'MID_PAC', 'BS_SAP', 'BS_LAC', 'SG_CLY', 'EXTRACTION_BLANK', 'LIBRARY_BLANK')

data.absolute$Sediment <- factor(data.absolute$Var2)
levels(data.absolute$Sediment) <- sedlist
data.absolute$Var1 <- factor(data.absolute$Var1, levels=c('157 bp to 177 bp', '177 bp to 202 bp', '202 bp to 227 bp', '227 bp to 277 bp', '277 bp to 327 bp', '327 bp to 427 bp'))

ggplot(data.absolute, aes(x=factor(Var2, level = level_order), y=value, fill=Var1)) + 
  geom_bar(position='stack', stat='identity', col="Black") +
  theme_bw() + 
  theme(legend.position = "right",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size=16, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(face="bold", size = 18, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=18),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=18),
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=14),
        strip.text.x = element_text(size = 14)) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,65)) +
  scale_fill_manual(values=rev(cols),guide = guide_legend(reverse = FALSE),limits = force) +
  ylab("Molarity [nM]") + xlab("Sample") +
  guides(fill=guide_legend(title="Fragment length")) +
  facet_grid(. ~ factor(Sediment, level = sedlevel),  scales = "free", space = "free")

#########################################################################################

############################ PLOTTING THE RATIO #########################################

ratio <- read.table("D:/All_the_data/Current_R_projects/EXTRACTIONS_ANALYSIS/FA_molarity_ratios_25-157-500.txt", sep="\t", header=T, stringsAsFactors = T, check.names = FALSE)

level_order <- c('400', '405', '412','417','401','406','413','418','402','407','414','419','403','408','415','420','404','409','416','421','410','411','422','423','Lib_blank_1(1-12)','Lib_blank_2(1-12)','Lib_blank_1(13-24)','Lib_blank_2(13-24)','Ladder') 

ggplot(ratio, aes(x=factor(name, level = level_order), y=ratio, shape=Extraction, fill=Buffer)) +
  geom_point(size = 2.5) +
  geom_linerange(aes(x=factor(name, level = level_order), ymax=ratio, ymin=0), position = "identity") +
  scale_shape_manual(values = c(23, 21, 22)) +
  scale_fill_manual(
    values = c("black", "blue", "red"),
    guide = guide_legend(override.aes = list(shape = 21))
  ) + theme_bw()

## you can also facet them per sediment type

sedlevel <- c('N-PAC', 'MID-PAC', 'BS-SAP', 'BS-LAC', 'SG-CLY', 'Extraction Blank', 'Library Blank')
level_order <- c('400', '412', '405','417','401','413','406','418','402','414','407','419','403','415','408','420','404','416','409','421','410','411','422','423','Lib_blank_1(1-12)','Lib_blank_2(1-12)','Lib_blank_1(13-24)','Lib_blank_2(13-24)') 

ggplot(ratio, aes(x=factor(name, level = level_order), y=ratio, shape=Extraction, fill=Buffer, group=name)) +
  geom_point(size = 3) + ylim(0,4) +
  geom_linerange(aes(x=factor(name, level = level_order), ymax=ratio, ymin=0, group=name), position = "identity") +
  scale_shape_manual(values = c(21, 22, 23), breaks=c('Magnetic Beads', 'Precipitation', 'Blank') ) +
  scale_fill_manual(
    values = c("red", "blue", "black"),
    breaks=c('Rohland', 'Direito', 'Blank'),
    guide = guide_legend(override.aes = list(shape = 21))) +
  theme(axis.text.x = element_text(size=12, angle = 90, vjust=0.4, hjust=1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(face="bold", size = 12, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=16),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=16),
        panel.background = element_rect(fill = "gray95"),
        panel.border = element_rect(colour = "black", fill = "NA"),
        legend.position = "right",
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.size = unit(0.5, "cm"),
        plot.title = element_text(size=16),
        strip.background = element_blank(), strip.text = element_text(size = 14, margin = margin(2,0,2,0))) +
  labs(title="Fragment molarity ratio", x ="Sample", y = "Molarity ratio [157-500 bp] / [25-157 bp]") +
  facet_grid(. ~ factor(Sediment, level = sedlevel), scales = "free", space = "free")
