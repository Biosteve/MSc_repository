# setwd("D:/Dropbox/DUTh/Thesis/SF") #laptop
setwd("C:/Users/User/Dropbox/DUTh/Thesis/SF") #PC

library(dplyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpval)
library(nortest)
library(Homo.sapiens)
library(GenomicRanges)
library(readr)
library(grid)
library(ComplexHeatmap)
library(UpSetR)
library(extrafont)
library(writexl)
library(RColorBrewer)

#===============================================================================

# remove the data that appear more than 4 times in Gene column in Average_0,1,7
#Maybe 3 times should be the ideal? 
# try with 2 times --> keep only the single ones

filtered_Average_0 <- Average_0 %>%
        group_by(Gene) %>%
        filter(n() < 2) %>%
        ungroup()

filtered_Average_1 <- Average_1 %>%
        group_by(Gene) %>%
        filter(n() < 2) %>%
        ungroup()

filtered_Average_7 <- Average_7 %>%
        group_by(Gene) %>%
        filter(n() < 2) %>%
        ungroup()


# UpSet Plots (similar to venn diagramms or Euler diagramms)

upset_0_low <- subset(filtered_Average_0, chr %in%
                              c('chr5', 'chr10', 'chr11', 'chr12', 'chr13',
                                'chr14', 'chr15', 'chr16', 'chr18', 'chr20', 'chr21'),
                      select = c(Gene))

upset_1_low <- subset(filtered_Average_1, chr %in% 
                              c('chr10', 'chr12', 'chr14', 'chr18', 'chr21'),
                      select = c(Gene))

upset_7_low <- subset(filtered_Average_7, chr %in% 
                              c('chr10', 'chr13', 'chr21'),
                      select = c(Gene))

list_low <- list(
        "Set_0_low" = unlist(upset_0_low),
        "Set_1_low" = unlist(upset_1_low),
        "Set_7_low" = unlist(upset_7_low)
)

Upset_low <- fromList(list_low)

upset(Upset_low,
      sets = c("Set_0_low", "Set_1_low", "Set_7_low"),
      nsets = 3, 
      nintersects = NA,
      order.by = "freq")

#=============================================================================
# 
# upset_0_mid <- subset(filtered_Average_0, chr %in%
#                               c('chr5', 'chr10', 'chr14', 'chr15', 'chr16'),
#                       select = c(Gene))
# 
# upset_1_mid <- subset(filtered_Average_1, chr %in% 
#                               c('chr1', 'chr2', 'chr22'),
#                       select = c(Gene))
# 
# upset_7_mid <- subset(filtered_Average_7, chr %in% 
#                               c('chr1', 'chr2', 'chr3', 'chr6', 'chr7', 'chr9',
#                                 'chr15', 'chr17', 'chr20', 'chr22'),
#                       select = c(Gene))
# 
# list_mid <- list(
#         "Set_0_mid" = unlist(upset_0_mid),
#         "Set_1_mid" = unlist(upset_1_mid),
#         "Set_7_mid" = unlist(upset_7_mid)
# )
# 
# Upset_mid <- fromList(list_mid)
# 
# upset(Upset_mid,
#       sets = c("Set_0_mid", "Set_1_mid", "Set_7_mid"),
#       nsets = 3, 
#       nintersects = NA,
#       order.by = "freq")

#=============================================================================

upset_0_high <- subset(filtered_Average_0, chr %in%
                              c('chr1', 'chr2', 'chr3', 'chr6', 'chr7', 'chr9',
                                'chr17', 'chr19', 'chr22'),
                      select = c(Gene))

upset_1_high <- subset(filtered_Average_1, chr %in% 
                              c('chr1', 'chr2', 'chr3', 'chr5', 'chr6', 'chr7',
                                'chr9', 'chr11','chr13', 'chr15', 'chr16',
                                'chr17', 'chr19', 'chr20', 'chr22'),
                      select = c(Gene))

upset_7_high <- subset(filtered_Average_7, chr %in% 
                              c('chr1', 'chr2', 'chr3', 'chr5',
                                'chr6', 'chr7', 'chr9', 'chr11', 'chr12',
                                'chr14', 'chr15','chr16', 'chr17', 'chr18',
                                'chr19', 'chr20','chr22'),
                      select = c(Gene))

list_high <- list(
        "Set_0_high" = unlist(upset_0_high),
        "Set_1_high" = unlist(upset_1_high),
        "Set_7_high" = unlist(upset_7_high)
)

Upset_high <- fromList(list_high)

upset(Upset_high,
      sets = c("Set_0_high", "Set_1_high", "Set_7_high"),
      nsets = 3, 
      nintersects = NA,
      order.by = "freq")

# ============================================================================

list_all <- list(
        "Set_0_low" = unlist(upset_0_low),
        "Set_1_low" = unlist(upset_1_low),
        "Set_7_low" = unlist(upset_7_low),
        "Set_0_high" = unlist(upset_0_high),
        "Set_1_high" = unlist(upset_1_high),
        "Set_7_high" = unlist(upset_7_high)
)

Upset_all <- fromList(list_all)

upset(Upset_all,
      sets = c("Set_7_high", "Set_7_low", "Set_1_high",
               "Set_1_low", "Set_0_high", "Set_0_low"),
      nsets = 6,
      nintersects = NA,
      keep.order = TRUE,
      show.numbers = 'Yes',
      scale.intersections = 'identity', 
      scale.sets = 'identity')

#==============================================================================
# 
# #Only the 0 & 7 day timepoints // best works with single copies of each gene!
# 
# list_0_7 <- list(
#         "Set_0_low" = unlist(upset_0_low),
#         "Set_7_low" = unlist(upset_7_low),
#         "Set_0_mid" = unlist(upset_0_mid),
#         "Set_7_mid" = unlist(upset_7_mid),
#         "Set_0_high" = unlist(upset_0_high),
#         "Set_7_high" = unlist(upset_7_high)
# )
# 
# Upset_0_7 <- fromList(list_0_7)
# 
# upset(Upset_0_7,
#       sets = c("Set_7_high", "Set_7_mid", "Set_7_low","Set_0_high", "Set_0_mid", "Set_0_low"),
#       nsets = 6,
#       nintersects = NA,
#       keep.order = TRUE,
#       show.numbers = 'Yes',
#       scale.intersections = 'identity', 
#       scale.sets = 'identity')


