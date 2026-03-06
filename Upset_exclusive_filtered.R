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

#==============================================================================

list_all <- list(
        "Set_0_low" = unlist(upset_0_low),
        "Set_1_low" = unlist(upset_1_low),
        "Set_7_low" = unlist(upset_7_low),
        "Set_0_high" = unlist(upset_0_high),
        "Set_1_high" = unlist(upset_1_high),
        "Set_7_high" = unlist(upset_7_high)
)

names(list_all) <- c("C (0 dpi)", "C (1 dpi)", "C (7 dpi)", 
                     "P (0 dpi)", "P (1 dpi)", "P (7 dpi)")


Upset_all <- fromList(list_all)

par(family = "serif")

tiff(filename = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/high_res_Upset_Plot.tiff",
     width = 12, height = 8, units = "in",
     res = 900, compression = "lzw", family = "serif")

upset(Upset_all,
      sets = c("P (7 dpi)", "C (7 dpi)", 
               "P (1 dpi)", "C (1 dpi)", 
               "P (0 dpi)", "C (0 dpi)"),
      nsets = 9,
      nintersects = NA,
      keep.order = TRUE,
      show.numbers = 'Yes',
      scale.intersections = 'identity', 
      scale.sets = 'identity',
      text.scale = c(2, 1.5, 2, 1.5, 2, 1.5), 
      point.size = 2.5,
      line.size = 0.8,
      mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "Number of Shared Genes", 
      sets.x.label = "Total Genes per Group",
      main.bar.color = "gray20", 
      sets.bar.color = "gray40",
      matrix.color = "black",
      shade.color = "gray90",
      )

dev.off()

all_genes <- unique(unlist(list_all))

all_genes <-as.data.frame(all_genes)
all_genes <- all_genes %>% mutate(ID_number = row_number()) %>%
        relocate(ID_number, all_genes)
all_genes$ID_number <- as.character(all_genes$ID_number)

gene_name_list <- data.frame(gene_name = all_genes, stringsAsFactors = FALSE)

#==============================================================================

# Exclusive 1304 genes among 3 groups

exclusive1304 <- Upset_all[
        Upset_all$"C (0 dpi)" == 1 & 
                Upset_all$"P (0 dpi)" == 0 &
                Upset_all$"C (1 dpi)" == 0 &
                Upset_all$"P (1 dpi)" == 1 &
                Upset_all$"C (7 dpi)" == 1 &
                Upset_all$"P (7 dpi)" == 0, 
]

gene_ID_1304 <- row.names(exclusive1304)
gene_1304 <- data.frame(gene_ID_1304)
colnames(gene_1304)[colnames(gene_1304) == "gene_ID_1304"] <- 'ID_number'
gene_1304 <- gene_1304 %>% left_join(all_genes, by='ID_number')

colnames(gene_1304)[colnames(gene_1304) == 'all_genes'] <- 'Gene'
gene_1304_list_0C_1P_7C <- merge(x=gene_1304, y=filtered_Average_0, by = 'Gene')
gene_1304_list_sum <- gene_1304_list_0C_1P_7C %>% count(Gene)

#==============================================================================

# Exclusive 2965 genes among 3 groups

exclusive2965 <- Upset_all[
        Upset_all$"C (0 dpi)" == 1 & 
                Upset_all$"P (0 dpi)" == 0 &
                Upset_all$"C (1 dpi)" == 1 &
                Upset_all$"P (1 dpi)" == 0 &
                Upset_all$"C (7 dpi)" == 1 &
                Upset_all$"P (7 dpi)" == 0, 
]

gene_ID_2965 <- row.names(exclusive2965)
gene_2965 <- data.frame(gene_ID_2965)
colnames(gene_2965)[colnames(gene_2965) == "gene_ID_2965"] <- 'ID_number'
gene_2965 <- gene_2965 %>% left_join(all_genes, by='ID_number')

colnames(gene_2965)[colnames(gene_2965) == 'all_genes'] <- 'Gene'
gene_2965_list_0C_1C_7C <- merge(x=gene_2965, y=filtered_Average_0, by = 'Gene')
gene_2965_list_sum <- gene_2965_list_0C_1C_7C %>% count(Gene)

#==============================================================================

# Exclusive 6193 genes among 3 groups

exclusive6193 <- Upset_all[
        Upset_all$"C (0 dpi)" == 1 & 
                Upset_all$"P (0 dpi)" == 0 &
                Upset_all$"C (1 dpi)" == 1 &
                Upset_all$"P (1 dpi)" == 0 &
                Upset_all$"C (7 dpi)" == 0 &
                Upset_all$"P (7 dpi)" == 1, 
]

gene_ID_6193 <- row.names(exclusive6193)
gene_6193 <- data.frame(gene_ID_6193)
colnames(gene_6193)[colnames(gene_6193) == "gene_ID_6193"] <- 'ID_number'
gene_6193 <- gene_6193 %>% left_join(all_genes, by='ID_number')

colnames(gene_6193)[colnames(gene_6193) == 'all_genes'] <- 'Gene'
gene_6193_list_0C_1C_7P <- merge(x=gene_6193, y=filtered_Average_0, by = 'Gene')
gene_6193_list_sum <- gene_6193_list_0C_1C_7P %>% count(Gene)


#==============================================================================

# Exclusive 11935 genes among 3 groups

exclusive11935 <- Upset_all[
        Upset_all$"C (0 dpi)" == 1 & 
                Upset_all$"P (0 dpi)" == 0 &
                Upset_all$"C (1 dpi)" == 0 &
                Upset_all$"P (1 dpi)" == 1 &
                Upset_all$"C (7 dpi)" == 0 &
                Upset_all$"P (7 dpi)" == 1, 
]

gene_ID_11935 <- row.names(exclusive11935)
gene_11935 <- data.frame(gene_ID_11935)
colnames(gene_11935)[colnames(gene_11935) == "gene_ID_11935"] <- 'ID_number'
gene_11935 <- gene_11935 %>% left_join(all_genes, by='ID_number')

colnames(gene_11935)[colnames(gene_11935) == 'all_genes'] <- 'Gene'
gene_11935_list_0C_1P_7P <- merge(x=gene_11935, y=filtered_Average_0, by = 'Gene')
gene_11935_list_sum <- gene_11935_list_0C_1P_7P %>% count(Gene)

#==============================================================================

# Exclusive 26999 genes among 3 groups

exclusive26999 <- Upset_all[
        Upset_all$"C (0 dpi)" == 0 & 
                Upset_all$"P (0 dpi)" == 1 &
                Upset_all$"C (1 dpi)" == 0 &
                Upset_all$"P (1 dpi)" == 1 &
                Upset_all$"C (7 dpi)" == 0 &
                Upset_all$"P (7 dpi)" == 1, 
]

gene_ID_26999 <- row.names(exclusive26999)
gene_26999 <- data.frame(gene_ID_26999)
colnames(gene_26999)[colnames(gene_26999) == "gene_ID_26999"] <- 'ID_number'
gene_26999 <- gene_26999 %>% left_join(all_genes, by='ID_number')

colnames(gene_26999)[colnames(gene_26999) == 'all_genes'] <- 'Gene'
gene_26999_list_0P_1P_7P <- merge(x=gene_26999, y=filtered_Average_0, by = 'Gene')
gene_26999_list_sum <- gene_26999_list_0P_1P_7P %>% count(Gene)


