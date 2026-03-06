#setwd("D:/Dropbox/DUTh/Thesis/SF") #laptop
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

data_desq2<-read.table("all_info_log2_normalised_counts_deseq2_scores_and_coordinates.txt")

colnames(data_desq2)[8] <- "Gene"

data_desq2$mean_Ctrl<-rowMeans(data_desq2[, c("Ctrl_1", "Ctrl_2")], na.rm = TRUE)
data_desq2$mean_H24<-rowMeans(data_desq2[, c("H24_1", "H24_2")], na.rm = TRUE)
data_desq2$mean_H168<-rowMeans(data_desq2[, c("H168_1", "H168_2")], na.rm = TRUE)


# Bring averages from Data_desq2 to the gene_1304 based on Genes

gene_1304_me <- merge(gene_1304, data_desq2, by='Gene') # n=272
gene_2965_me <- merge(gene_2965, data_desq2, by='Gene') # n=738 
gene_6193_me <- merge(gene_6193, data_desq2, by='Gene') # n=1595
gene_11935_me <- merge(gene_11935, data_desq2, by='Gene') # n=3223
gene_26999_me <- merge(gene_26999, data_desq2, by='Gene') # n=8084

# Create the box plots for each cluster

my_comparisons <- list(c('mean_Ctrl', 'mean_H168'), c('mean_Ctrl', 'mean_H24'),
                       c('mean_H24', 'mean_H168'))

#=========1304================================================================

data_1304 <- gene_1304_me %>%
        pivot_longer(cols=c(mean_Ctrl, mean_H24, mean_H168), 
                     names_to = 'variable',values_to = 'value')

final_plot_1304_cluster <- ggplot(data_1304, aes(x = variable, y = value, fill = variable)) +
        scale_y_continuous(limits = c(0,20)) +
        geom_boxplot(alpha=1) +
        stat_compare_means(comparisons= my_comparisons,
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="1304 cluster (n=272)") + 
        labs(x= "Time points", y="Average distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggsave("high_res_final_ggplot_1304_cluster.tiff",
       plot = final_plot_1304_cluster,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


#=========2965================================================================

data_2965 <- gene_2965_me %>%
        pivot_longer(cols=c(mean_Ctrl, mean_H24, mean_H168), 
                     names_to = 'variable',values_to = 'value')

final_plot_2965_cluster <- ggplot(data_2965, aes(x = variable, y = value, fill = variable)) +
        scale_y_continuous(limits = c(0,20)) +
        geom_boxplot(alpha=1) +
        stat_compare_means(comparisons= my_comparisons,
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="2965 cluster (n=738)") + 
        labs(x= "Time points", y="Average distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggsave("high_res_final_ggplot_2965_cluster.tiff",
       plot = final_plot_2965_cluster,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


#=========6193================================================================

data_6193 <- gene_6193_me %>%
        pivot_longer(cols=c(mean_Ctrl, mean_H24, mean_H168), 
                     names_to = 'variable',values_to = 'value')

final_plot_6193_cluster <- ggplot(data_6193, aes(x = variable, y = value, fill = variable)) +
        scale_y_continuous(limits = c(0,25)) +
        geom_boxplot(alpha=1) +
        stat_compare_means(comparisons= my_comparisons,
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="6193 cluster (n=1595)") + 
        labs(x= "Time points", y="Average distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggsave("high_res_final_ggplot_6193_cluster.tiff",
       plot = final_plot_6193_cluster,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')



#=========11935================================================================

data_11935 <- gene_11935_me %>%
        pivot_longer(cols=c(mean_Ctrl, mean_H24, mean_H168), 
                     names_to = 'variable',values_to = 'value')

final_plot_11935_cluster <- ggplot(data_11935, aes(x = variable, y = value, fill = variable)) +
        scale_y_continuous(limits = c(0,25)) +
        geom_boxplot(alpha=1) +
        stat_compare_means(comparisons= my_comparisons,
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="11935 cluster (n=3223)") + 
        labs(x= "Time points", y="Average distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggsave("high_res_final_ggplot_11935_cluster.tiff",
       plot = final_plot_11935_cluster,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')



#=========26999================================================================

data_26999 <- gene_26999_me %>%
        pivot_longer(cols=c(mean_Ctrl, mean_H24, mean_H168), 
                     names_to = 'variable',values_to = 'value')

final_plot_26999_cluster <- ggplot(data_26999, aes(x = variable, y = value, fill = variable)) +
        scale_y_continuous(limits = c(0,30)) +
        geom_boxplot(alpha=1) +
        stat_compare_means(comparisons= my_comparisons,
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="26999 cluster (n=8084)") + 
        labs(x= "Time points", y="Average distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggsave("high_res_final_ggplot_26999_cluster.tiff",
       plot = final_plot_26999_cluster,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

