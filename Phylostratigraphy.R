setwd("D:/Dropbox/DUTh/Thesis/SF") #laptop
#setwd("C:/Users/User/Dropbox/DUTh/Thesis/SF") #PC

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

##Phylostratigraphy##

Phylo_data<-read.table("HSapiens_Phylostratigraphy_Genes_Annotated_hg19coords.bed")

summary(Phylo_data)
test_phylo<-Phylo_data[1:9]
sum(is.na(test_phylo))

unique(Phylo_data$V9) # Gives 16 species
colnames(Phylo_data)=c("chr", "startPOS", "endPOS", "Gene", "+/-", "Ensembl", "NCBI_Gene", "Spec_No", "Species")

## Rename the columns // subset based on species // Find the genes in each one 

LUCA_set <- Phylo_data[Phylo_data$"Species" == 'LUCA',]
Sapiens_set <- Phylo_data[Phylo_data$"Species" == 'HomoSapiens',]
Eutheria_set <- Phylo_data[Phylo_data$"Species" == 'Eutheria',]
Eukaryota_set <- Phylo_data[Phylo_data$"Species" == 'Eukaryota',]
Metazoa_set <- Phylo_data[Phylo_data$"Species" == 'Metazoa',]
Euteleostomi_set <- Phylo_data[Phylo_data$"Species" == 'Euteleostomi',]
Theria_set <- Phylo_data[Phylo_data$"Species" == 'Theria',]
Eumetazoa_set <- Phylo_data[Phylo_data$"Species" == 'Eumetazoa',]
Homininae_set <- Phylo_data[Phylo_data$"Species" == 'Homininae',]
Amniota_set <- Phylo_data[Phylo_data$"Species" == 'Amniota',]
Mammalia_set <- Phylo_data[Phylo_data$"Species" == 'Mammalia',]
Chordata_set <- Phylo_data[Phylo_data$"Species" == 'Chordata',]
Bilateria_set <- Phylo_data[Phylo_data$"Species" == 'Bilateria',]
Euarchontoglires_set <- Phylo_data[Phylo_data$"Species" == 'Euarchontoglires',]
Opisthokonta_set <- Phylo_data[Phylo_data$"Species" == 'Opisthokonta',]
Catarrhini_set <- Phylo_data[Phylo_data$"Species" == 'Catarrhini',]

## Function for automated subset
# 
# subset_by_column <- function (data, column, value){
#         data[data[[column]] == value,]
# }
# 
# Metazoa_set=subset_by_column(Phylo_data, 'V9', 'Metazoa')

##Use HomoSapiens, with small number as a template for the rest species##


## Must become a function to reduce the derivatives and
## automate the final result

## For Homo Sapiens ##

# Sapiens_set2 <- Sapiens_set %>% inner_join(Average_0, by = "Gene")
# Sapiens_set2 <- Sapiens_set2[, -c(10:19)]
# colnames(Sapiens_set2)[colnames(Sapiens_set2) == "average"] <- "Average_0"
# 
# Sapiens_set3 <- Sapiens_set2 %>% inner_join(Average_1, by = "Gene")
# Sapiens_set3 <- Sapiens_set3[, -c(11:20)]
# colnames(Sapiens_set3)[colnames(Sapiens_set3) == "average"] <- "Average_1"
# 
# Sapiens_set4 <- Sapiens_set3 %>% inner_join(Average_7, by = "Gene")
# Sapiens_set4 <- Sapiens_set4[, -c(12:21)]
# colnames(Sapiens_set4)[colnames(Sapiens_set4) == "average"] <- "Average_7"
# 
# Sapiens_set4 <- na.omit(Sapiens_set4) #removes rows with NaNs
# 
# graphing_SS4_0<- ggplot(Sapiens_set4, aes(x=Gene, y=Average_0, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_SS4_1<- ggplot(Sapiens_set4, aes(x=Gene, y=Average_1, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_SS4_7<- ggplot(Sapiens_set4, aes(x=Gene, y=Average_7, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# 
# ggarrange(graphing_SS4_0, graphing_SS4_1, graphing_SS4_7, ncol =3 , nrow = 1)


# Sapiens_tibble <- as_tibble(Sapiens)
# Sapiens_tibble %>% group_by(chr.x) %>% summarise(n = n())


## For LUCA ##
# 
# LUCA_set2 <- LUCA_set %>% inner_join(Average_0, by = "Gene")
# LUCA_set2 <- LUCA_set2[, -c(10:19)]
# colnames(LUCA_set2)[colnames(LUCA_set2) == "average"] <- "Average_0"
# 
# LUCA_set3 <- LUCA_set2 %>% inner_join(Average_1, by = "Gene")
# LUCA_set3 <- LUCA_set3[, -c(11:20)]
# colnames(LUCA_set3)[colnames(LUCA_set3) == "average"] <- "Average_1"
# 
# LUCA_set4 <- LUCA_set3 %>% inner_join(Average_7, by = "Gene")
# LUCA_set4 <- LUCA_set4[, -c(12:21)]
# colnames(LUCA_set4)[colnames(LUCA_set4) == "average"] <- "Average_7"
# 
# LUCA_set4 <- na.omit(LUCA_set4) #removes rows with NaNs

# 
# graphing_LS4_0<- ggplot(LUCA_set4, aes(x=Gene, y=Average_0, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_LS4_1<- ggplot(LUCA_set4, aes(x=Gene, y=Average_1, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_LS4_7<- ggplot(LUCA_set4, aes(x=Gene, y=Average_7, color=Gene)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Testing") + labs(x= "Gene", y="Mean Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# 
# ggarrange(graphing_LS4_0, graphing_LS4_1, graphing_LS4_7, ncol =3 , nrow = 1)
# 
# LUCA_tibble <- as_tibble(LUCA)
# LUCA_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)


## For Mammalia ##
# 
# Mammalia_set2 <- Mammalia_set %>% inner_join(Average_0, by = "Gene")
# Mammalia_set2 <- Mammalia_set2[, -c(10:19)]
# colnames(Mammalia_set2)[colnames(Mammalia_set2) == "average"] <- "Average_0"
# 
# Mammalia_set3 <- Mammalia_set2 %>% inner_join(Average_1, by = "Gene")
# Mammalia_set3 <- Mammalia_set3[, -c(11:20)]
# colnames(Mammalia_set3)[colnames(Mammalia_set3) == "average"] <- "Average_1"
# 
# Mammalia_set4 <- Mammalia_set3 %>% inner_join(Average_7, by = "Gene")
# Mammalia_set4 <- Mammalia_set4[, -c(12:21)]
# colnames(Mammalia_set4)[colnames(Mammalia_set4) == "average"] <- "Average_7"
# 
# Mammalia_set4 <- na.omit(Mammalia_set4) #removes rows with NaNs
# 
# 
# Mammalia_tibble <- as_tibble(Mammalia)
# Mammalia_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)


## Functioning & Sub-setting 

subsetting <- function (x){
        x.1 <- x %>% inner_join(Average_0, by = "Gene")
                x.1 <- x.1[, -c(10:19)]
                colnames(x.1)[colnames(x.1) == "average"] <- "Average_0"
        x.2 <- x.1 %>% inner_join(Average_1, by = "Gene")
                x.2 <- x.2[, -c(11:20)]
                colnames(x.2)[colnames(x.2) == "average"] <- "Average_1"
        x.3 <- x.2 %>% inner_join(Average_7, by = "Gene")
                x.3 <- x.3[, -c(12:21)]
                colnames(x.3)[colnames(x.3) == "average"] <- "Average_7"
        x.4 <- x.3 %>% mutate(Total_average=rowMeans(x.3[, c('Average_0', 'Average_1', 'Average_7' )]))
        x.5 <- na.omit(x.4)
        return(x.5)
}

##============================================================================##

Sapiens <- subsetting(Sapiens_set)
Sapiens_tibble <- as_tibble(Sapiens)
Sapiens_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Mammalia <- subsetting(Mammalia_set)
Mammalia_tibble <- as_tibble(Mammalia)
Mammalia_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

LUCA <- subsetting(LUCA_set)
LUCA_tibble <- as_tibble(LUCA)
LUCA_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Eutheria <- subsetting(Eutheria_set)
Eutheria_tibble <- as_tibble(Eutheria)
Eutheria_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Eukaryota <- subsetting(Eukaryota_set)
Eukaryota_tibble <- as_tibble(Eukaryota)
Eukaryota_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Metazoa <- subsetting(Metazoa_set)
Metazoa_tibble <- as_tibble(Metazoa)
Metazoa_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Euteleostomi <- subsetting(Euteleostomi_set)
Euteleostomi_tibble <- as_tibble(Euteleostomi)
Euteleostomi_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Theria <- subsetting(Theria_set)
Theria_tibble <- as_tibble(Theria)
Theria_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Eumetazoa <- subsetting(Eumetazoa_set)
Eumetazoa_tibble <- as_tibble(Eumetazoa)
Eumetazoa_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Homininae <- subsetting(Homininae_set)
Homininae_tibble <- as_tibble(Homininae)
Homininae_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Amniota <- subsetting(Amniota_set)
Amniota_tibble <- as_tibble(Amniota)
Amniota_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Chordata <- subsetting(Chordata_set)
Chordata_tibble <- as_tibble(Chordata)
Chordata_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Bilateria <- subsetting(Bilateria_set)
Bilateria_tibble <- as_tibble(Bilateria)
Bilateria_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Euarchontoglires <- subsetting(Euarchontoglires_set)
Euarchontoglires_tibble <- as_tibble(Euarchontoglires)
Euarchontoglires_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Opisthokonta <- subsetting(Opisthokonta_set)
Opisthokonta_tibble <- as_tibble(Opisthokonta)
Opisthokonta_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)

Catarrhini <- subsetting(Catarrhini_set)
Catarrhini_tibble <- as_tibble(Catarrhini)
Catarrhini_tibble %>% group_by(chr.x) %>% summarise(n = n()) %>% print(n=Inf)


##============================================================================##

## Check for chrs in Sapiens

Sapiens_sum <- Sapiens %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in LUCA

LUCA_sum <- LUCA %>% group_by(chr.x) %>%
        filter(chr.x != "chrUn_gl000212") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

#write_xlsx(LUCA_sum, 'LUCA_sum_data.xlsx')


## Check for chrs in Mammalia

Mammalia_sum <- Mammalia %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

#write_xlsx(Mammalia_sum, 'Mammalia_sum_data.xlsx')


## Check for chrs in Eutheria

Eutheria_sum <- Eutheria %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Eukaryota

Eukaryota_sum <- Eukaryota %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Metazoa

Metazoa_sum <- Metazoa %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Euteleostomi

Euteleostomi_sum <- Euteleostomi %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Theria

Theria_sum <- Theria %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Eumetazoa

Eumetazoa_sum <- Eumetazoa %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Homininae

Homininae_sum <- Homininae %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Amniota

Amniota_sum <- Amniota %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Chordata

Chordata_sum <- Chordata %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Bilateria

Bilateria_sum <- Bilateria %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Euarchontoglires

Euarchontoglires_sum <- Euarchontoglires %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Opisthokonta

Opisthokonta_sum <- Opisthokonta %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

## Check for chrs in Catarrhini

Catarrhini_sum <- Catarrhini %>% group_by(chr.x) %>%
        filter(chr.x != "NA") %>%
        summarise(across(c(Average_0, Average_1, Average_7, Total_average), mean), .groups = 'drop') %>%
        as.data.frame()

##============================================================================##

correct_order <- c("chr1", "chr2", "chr3", "chr5", "chr6", "chr7", "chr9", "chr10", 
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                   "chr18", "chr19", "chr20", "chr21", "chr22")


# LUCA graphs

LUCA_sum$chr.x <- factor(LUCA_sum$chr.x)

LUCA$chr.x <- factor(LUCA$chr.x, levels = correct_order)

LUCA_graph_0 <-ggplot(subset(LUCA, !is.na(chr.x)), aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_LUCA_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = LUCA_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

LUCA_graph_1 <-ggplot(subset(LUCA, !is.na(chr.x)), aes(x = chr.x, y = Average_1))  + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_LUCA_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = LUCA_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

LUCA_graph_7 <-ggplot(subset(LUCA, !is.na(chr.x)), aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_LUCA_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = LUCA_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(LUCA_graph_0, LUCA_graph_1, LUCA_graph_7, ncol =1 , nrow = 3)

combined_plot_LUCA <- ggarrange(LUCA_graph_0, LUCA_graph_1, LUCA_graph_7, ncol =1 , nrow = 3)

final_plot_LUCA <- annotate_figure(
        combined_plot_LUCA,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif',
                                                 fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90,
                        gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_LUCA)

ggsave("high_res_final_plot_LUCA.tiff",
       plot = final_plot_LUCA,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Sapiens graphs

Sapiens_sum$chr.x <- factor(Sapiens_sum$chr.x)

Sapiens$chr.x <- factor(Sapiens$chr.x, levels = correct_order)

Sapiens_graph_0 <-ggplot(Sapiens, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Sapiens_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Sapiens_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Sapiens_graph_1 <-ggplot(Sapiens, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Sapiens_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Sapiens_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Sapiens_graph_7 <-ggplot(Sapiens, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Sapiens_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Sapiens_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Sapiens_graph_0, Sapiens_graph_1, Sapiens_graph_7, ncol =1 , nrow = 3)

combined_plot_Sapiens <- ggarrange(Sapiens_graph_0, Sapiens_graph_1, Sapiens_graph_7, ncol =1 , nrow = 3)

final_plot_Sapiens <- annotate_figure(
        combined_plot_Sapiens,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif',
                                                 fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90,
                        gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Sapiens)

ggsave("high_res_final_plot_Sapiens.tiff",
       plot = final_plot_Sapiens,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Mammalia graphs

Mammalia_sum$chr.x <- factor(Mammalia_sum$chr.x)

Mammalia$chr.x <- factor(Mammalia$chr.x, levels = correct_order)


Mammalia_graph_0 <-ggplot(Mammalia, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Mammalia_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)


Mammalia_graph_1 <-ggplot(Mammalia, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Mammalia_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Mammalia_graph_7 <-ggplot(Mammalia, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Mammalia_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)


combined_plot_Mammalia <- ggarrange(Mammalia_graph_0, Mammalia_graph_1, Mammalia_graph_7, ncol =1 , nrow = 3)

final_plot_Mammalia <- annotate_figure(
        combined_plot_Mammalia,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif',
                                                 fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90,
                        gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Mammalia)

ggsave("high_res_final_plot_Mammalia.tiff",
       plot = final_plot_Mammalia,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')



#gia overall graph x= "Chromosomes", y="Distance from Center"

# Eutheria graphs

Eutheria_sum$chr.x <- factor(Eutheria_sum$chr.x)

Eutheria$chr.x <- factor(Eutheria$chr.x, levels = correct_order)

Eutheria_graph_0 <-ggplot(Eutheria, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eutheria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)


Eutheria_graph_1 <-ggplot(Eutheria, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eutheria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Eutheria_graph_7 <-ggplot(Eutheria, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) +
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eutheria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Eutheria_graph_0, Eutheria_graph_1, Eutheria_graph_7, ncol =1 , nrow = 3)


combined_plot_Eutheria <- ggarrange(Eutheria_graph_0, Eutheria_graph_1, Eutheria_graph_7, ncol =1 , nrow = 3)


final_plot_Eutheria <- annotate_figure(
        combined_plot_Eutheria,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Eutheria)

ggsave("high_res_final_plot_Eutheria.tiff",
       plot = final_plot_Eutheria,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Eukaryota graphs

Eukaryota_sum$chr.x <- factor(Eukaryota_sum$chr.x)

Eukaryota$chr.x <- factor(Eukaryota$chr.x, levels = correct_order)

Eukaryota_graph_0 <-ggplot(Eukaryota, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eukaryota_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eukaryota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Eukaryota_graph_1 <-ggplot(Eukaryota, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eukaryota_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eukaryota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Eukaryota_graph_7 <-ggplot(Eukaryota, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eukaryota_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eukaryota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Eukaryota_graph_0, Eukaryota_graph_1, Eukaryota_graph_7, ncol =1 , nrow = 3)

combined_plot_Eukaryota <- ggarrange(Eukaryota_graph_0, Eukaryota_graph_1, Eukaryota_graph_7, ncol =1 , nrow = 3)


final_plot_Eukaryota <- annotate_figure(
        combined_plot_Eukaryota,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Eukaryota)

ggsave("high_res_final_plot_Eukaryota.tiff",
       plot = final_plot_Eukaryota,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Metazoa graphs

Metazoa_sum$chr.x <- factor(Metazoa_sum$chr.x)

Metazoa$chr.x <- factor(Metazoa$chr.x, levels = correct_order)


Metazoa_graph_0 <-ggplot(Metazoa, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Metazoa_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Metazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Metazoa_graph_1 <-ggplot(Metazoa, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Metazoa_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Metazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Metazoa_graph_7 <-ggplot(Metazoa, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Metazoa_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Metazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Metazoa_graph_0, Metazoa_graph_1, Metazoa_graph_7, ncol =1 , nrow = 3)

combined_plot_Metazoa <- ggarrange(Metazoa_graph_0, Metazoa_graph_1, Metazoa_graph_7, ncol =1 , nrow = 3)


final_plot_Metazoa <- annotate_figure(
        combined_plot_Metazoa,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Metazoa)

ggsave("high_res_final_plot_Metazoa.tiff",
       plot = final_plot_Metazoa,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Euteleostomi graphs

Euteleostomi_sum$chr.x <- factor(Euteleostomi_sum$chr.x)

Euteleostomi$chr.x <- factor(Euteleostomi$chr.x, levels=correct_order)


Euteleostomi_graph_0 <-ggplot(Euteleostomi, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euteleostomi_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euteleostomi_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Euteleostomi_graph_1 <-ggplot(Euteleostomi, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euteleostomi_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euteleostomi_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Euteleostomi_graph_7 <-ggplot(Euteleostomi, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euteleostomi_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euteleostomi_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Euteleostomi_graph_0, Euteleostomi_graph_1, Euteleostomi_graph_7, ncol =1 , nrow = 3)

combined_plot_Euteleostomi <- ggarrange(Euteleostomi_graph_0, Euteleostomi_graph_1, Euteleostomi_graph_7, ncol =1 , nrow = 3)


final_plot_Euteleostomi <- annotate_figure(
        combined_plot_Euteleostomi,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Euteleostomi)

ggsave("high_res_final_plot_Euteleostomi.tiff",
       plot = final_plot_Euteleostomi,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Theria graphs

Theria_sum$chr.x <- factor(Theria_sum$chr.x)

Theria$chr.x <- factor(Theria$chr.x, levels=correct_order)

Theria_graph_0 <-ggplot(Theria, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Theria_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Theria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Theria_graph_1 <-ggplot(Theria, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Theria_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Theria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Theria_graph_7 <-ggplot(Theria, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Theria_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Theria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Theria_graph_0, Theria_graph_1, Theria_graph_7, ncol =1 , nrow = 3)

combined_plot_Theria <- ggarrange(Theria_graph_0, Theria_graph_1, Theria_graph_7, ncol =1 , nrow = 3)


final_plot_Theria <- annotate_figure(
        combined_plot_Theria,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Theria)

ggsave("high_res_final_plot_Theria.tiff",
       plot = final_plot_Theria,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Eumetazoa graphs

Eumetazoa_sum$chr.x <- factor(Eumetazoa_sum$chr.x)

Eumetazoa$chr.x <- factor(Eumetazoa$chr.x, levels=correct_order)


Eumetazoa_graph_0 <-ggplot(Eumetazoa, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eumetazoa_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eumetazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Eumetazoa_graph_1 <-ggplot(Eumetazoa, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eumetazoa_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eumetazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Eumetazoa_graph_7 <-ggplot(Eumetazoa, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Eumetazoa_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Eumetazoa_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Eumetazoa_graph_0, Eumetazoa_graph_1, Eumetazoa_graph_7, ncol =1 , nrow = 3)


combined_plot_Eumetazoa <- ggarrange(Eumetazoa_graph_0, Eumetazoa_graph_1, Eumetazoa_graph_7, ncol =1 , nrow = 3)


final_plot_Eumetazoa <- annotate_figure(
        combined_plot_Eumetazoa,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Eumetazoa)

ggsave("high_res_final_plot_Eumetazoa.tiff",
       plot = final_plot_Eumetazoa,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Homininae graphs

Homininae_sum$chr.x <- factor(Homininae_sum$chr.x)

Homininae$chr.x <- factor(Homininae$chr.x, levels=correct_order)


Homininae_graph_0 <-ggplot(Homininae, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Homininae_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Homininae_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Homininae_graph_1 <-ggplot(Homininae, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Homininae_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Homininae_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Homininae_graph_7 <-ggplot(Homininae, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Homininae_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Homininae_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Homininae_graph_0, Homininae_graph_1, Homininae_graph_7, ncol =1 , nrow = 3)

combined_plot_Homininae <- ggarrange(Homininae_graph_0, Homininae_graph_1, Homininae_graph_7, ncol =1 , nrow = 3)

final_plot_Homininae <- annotate_figure(
        combined_plot_Homininae,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Homininae)

ggsave("high_res_final_plot_Homininae.tiff",
       plot = final_plot_Homininae,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Amniota graphs

Amniota_sum$chr.x <- factor(Amniota_sum$chr.x)

Amniota$chr.x <- factor(Amniota$chr.x, levels=correct_order)

Amniota_graph_0 <-ggplot(Amniota, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Amniota_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Amniota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Amniota_graph_1 <-ggplot(Amniota, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Amniota_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Amniota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Amniota_graph_7 <-ggplot(Amniota, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Amniota_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Amniota_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Amniota_graph_0, Amniota_graph_1, Amniota_graph_7, ncol =1 , nrow = 3)


combined_plot_Amniota <- ggarrange(Amniota_graph_0, Amniota_graph_1, Amniota_graph_7, ncol =1 , nrow = 3)

final_plot_Amniota <- annotate_figure(
        combined_plot_Amniota,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Amniota)

ggsave("high_res_final_plot_Amniota.tiff",
       plot = final_plot_Amniota,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Chordata graphs

Chordata_sum$chr.x <- factor(Chordata_sum$chr.x)

Chordata$chr.x <- factor(Chordata$chr.x, levels=correct_order)

Chordata_graph_0 <-ggplot(Chordata, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Chordata_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Chordata_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Chordata_graph_1 <-ggplot(Chordata, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Chordata_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Chordata_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Chordata_graph_7 <-ggplot(Chordata, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Chordata_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Chordata_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Chordata_graph_0, Chordata_graph_1, Chordata_graph_7, ncol =1 , nrow = 3)

combined_plot_Chordata <- ggarrange(Chordata_graph_0, Chordata_graph_1, Chordata_graph_7, ncol =1 , nrow = 3)

final_plot_Chordata <- annotate_figure(
        combined_plot_Chordata,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Chordata)

ggsave("high_res_final_plot_Chordata.tiff",
       plot = final_plot_Chordata,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Bilateria graphs

Bilateria_sum$chr.x <- factor(Bilateria_sum$chr.x)

Bilateria$chr.x <- factor(Bilateria$chr.x, levels=correct_order)


Bilateria_graph_0 <-ggplot(Bilateria, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Bilateria_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Bilateria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Bilateria_graph_1 <-ggplot(Bilateria, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Bilateria_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Bilateria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Bilateria_graph_7 <-ggplot(Bilateria, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Bilateria_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Bilateria_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Bilateria_graph_0, Bilateria_graph_1, Bilateria_graph_7, ncol =1 , nrow = 3)

combined_plot_Bilateria <- ggarrange(Bilateria_graph_0, Bilateria_graph_1, Bilateria_graph_7, ncol =1 , nrow = 3)

final_plot_Bilateria <- annotate_figure(
        combined_plot_Bilateria,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Bilateria)

ggsave("high_res_final_plot_Bilateria.tiff",
       plot = final_plot_Bilateria,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Euarchontoglires graphs

Euarchontoglires_sum$chr.x <- factor(Euarchontoglires_sum$chr.x)

Euarchontoglires$chr.x <- factor(Euarchontoglires$chr.x, levels=correct_order)

Euarchontoglires_graph_0 <-ggplot(Euarchontoglires, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euarchontoglires_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euarchontoglires_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Euarchontoglires_graph_1 <-ggplot(Euarchontoglires, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euarchontoglires_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euarchontoglires_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Euarchontoglires_graph_7 <-ggplot(Euarchontoglires, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Euarchontoglires_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Euarchontoglires_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Euarchontoglires_graph_0, Euarchontoglires_graph_1, Euarchontoglires_graph_7, ncol =1 , nrow = 3)

combined_plot_Euarchontoglires <- ggarrange(Euarchontoglires_graph_0, Euarchontoglires_graph_1, Euarchontoglires_graph_7, ncol =1 , nrow = 3)

final_plot_Euarchontoglires <- annotate_figure(
        combined_plot_Euarchontoglires,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Euarchontoglires)

ggsave("high_res_final_plot_Euarchontoglires.tiff",
       plot = final_plot_Euarchontoglires,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# Opisthokonta graphs

Opisthokonta_sum$chr.x <- factor(Opisthokonta_sum$chr.x)

Opisthokonta$chr.x <- factor(Opisthokonta$chr.x, levels=correct_order)

Opisthokonta_graph_0 <-ggplot(Opisthokonta, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Opisthokonta_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Opisthokonta_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Opisthokonta_graph_1 <-ggplot(Opisthokonta, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Opisthokonta_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Opisthokonta_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Opisthokonta_graph_7 <-ggplot(Opisthokonta, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Opisthokonta_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Opisthokonta_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Opisthokonta_graph_0, Opisthokonta_graph_1, Opisthokonta_graph_7, ncol =1 , nrow = 3)

combined_plot_Opisthokonta <- ggarrange(Opisthokonta_graph_0, Opisthokonta_graph_1, Opisthokonta_graph_7, ncol =1 , nrow = 3)

final_plot_Opisthokonta <- annotate_figure(
        combined_plot_Opisthokonta,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Opisthokonta)

ggsave("high_res_final_plot_Opisthokonta.tiff",
       plot = final_plot_Opisthokonta,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# Catarrhini graphs

Catarrhini_sum$chr.x <- factor(Catarrhini_sum$chr.x)

Catarrhini$chr.x <- factor(Catarrhini$chr.x, levels=correct_order)

Catarrhini_graph_0 <-ggplot(Catarrhini, aes(x = chr.x, y = Average_0)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Catarrhini_chr_Dist0") + 
        labs(title="0 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Catarrhini_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Catarrhini_graph_1 <-ggplot(Catarrhini, aes(x = chr.x, y = Average_1)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Catarrhini_chr_Dist1") + 
        labs(title="1 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Catarrhini_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

Catarrhini_graph_7 <-ggplot(Catarrhini, aes(x = chr.x, y = Average_7)) + 
        scale_x_discrete() +
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="Mean_Catarrhini_chr_Dist7") + 
        labs(title="7 dpi") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = 'black'),
              axis.title.x = element_blank(),
              axis.line.y = element_line(linewidth = 1, colour = 'black'),
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank()) + 
        geom_segment(data = Catarrhini_sum,
                     aes(x = as.numeric(chr.x) - 0.4, 
                         xend = as.numeric(chr.x) + 0.4, 
                         y = Total_average, yend = Total_average), 
                     inherit.aes = TRUE, colour = 'orange', size = 1)

ggarrange(Catarrhini_graph_0, Catarrhini_graph_1, Catarrhini_graph_7, ncol =1 , nrow = 3)

combined_plot_Catarrhini <- ggarrange(Catarrhini_graph_0, Catarrhini_graph_1, Catarrhini_graph_7, ncol =1 , nrow = 3)

final_plot_Catarrhini <- annotate_figure(
        combined_plot_Catarrhini,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distance from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Catarrhini)

ggsave("high_res_final_plot_Catarrhini.tiff",
       plot = final_plot_Catarrhini,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Phylo",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')

# first we create the data frame for the heatmap

Phylo_sum <- data.frame(chr.x = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                                  'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                  'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                  'chr21', 'chr22', 'chrX'))

# Data table for d0 dataset

Phylo_sum_0 <- Phylo_sum %>%
        left_join(dplyr::select(LUCA_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'LUCA_0') %>%
        left_join(dplyr::select(Mammalia_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Mammalia_0') %>%
        left_join(dplyr::select(Sapiens_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Sapiens_0') %>%
        left_join(dplyr::select(Eutheria_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Eutheria_0') %>%
        left_join(dplyr::select(Eukaryota_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Eukaryota_0') %>%
        left_join(dplyr::select(Metazoa_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Metazoa_0') %>%
        left_join(dplyr::select(Euteleostomi_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Euteleostomi_0') %>%
        left_join(dplyr::select(Theria_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Theria_0') %>%
        left_join(dplyr::select(Eumetazoa_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Eumetazoa_0') %>%
        left_join(dplyr::select(Homininae_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Homininae_0') %>%
        left_join(dplyr::select(Amniota_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Amniota_0') %>%
        left_join(dplyr::select(Chordata_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Chordata_sum_0') %>%
        left_join(dplyr::select(Bilateria_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Bilateria_0') %>%
        left_join(dplyr::select(Euarchontoglires_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Euarchontoglires_0') %>%
        left_join(dplyr::select(Opisthokonta_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Opisthokonta_0') %>%
        left_join(dplyr::select(Catarrhini_sum, chr.x, Average_0), by = 'chr.x') %>%
        rename('Average_0' = 'Catarrhini_0')

# Data table for d1 dataset

Phylo_sum_1 <- Phylo_sum %>%
        left_join(dplyr::select(LUCA_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'LUCA_1') %>%
        left_join(dplyr::select(Mammalia_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Mammalia_1') %>%
        left_join(dplyr::select(Sapiens_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Sapiens_1') %>%
        left_join(dplyr::select(Eutheria_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Eutheria_1') %>%
        left_join(dplyr::select(Eukaryota_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Eukaryota_1') %>%
        left_join(dplyr::select(Metazoa_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Metazoa_1') %>%
        left_join(dplyr::select(Euteleostomi_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Euteleostomi_1') %>%
        left_join(dplyr::select(Theria_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Theria_1') %>%
        left_join(dplyr::select(Eumetazoa_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Eumetazoa_1') %>%
        left_join(dplyr::select(Homininae_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Homininae_1') %>%
        left_join(dplyr::select(Amniota_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Amniota_1') %>%
        left_join(dplyr::select(Chordata_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Chordata_sum_1') %>%
        left_join(dplyr::select(Bilateria_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Bilateria_1') %>%
        left_join(dplyr::select(Euarchontoglires_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Euarchontoglires_1') %>%
        left_join(dplyr::select(Opisthokonta_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Opisthokonta_1') %>%
        left_join(dplyr::select(Catarrhini_sum, chr.x, Average_1), by = 'chr.x') %>%
        rename('Average_1' = 'Catarrhini_1')


# Data table for d7 dataset

Phylo_sum_7 <- Phylo_sum %>%
        left_join(dplyr::select(LUCA_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'LUCA_7') %>%
        left_join(dplyr::select(Mammalia_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Mammalia_7') %>%
        left_join(dplyr::select(Sapiens_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Sapiens_7') %>%
        left_join(dplyr::select(Eutheria_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Eutheria_7') %>%
        left_join(dplyr::select(Eukaryota_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Eukaryota_7') %>%
        left_join(dplyr::select(Metazoa_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Metazoa_7') %>%
        left_join(dplyr::select(Euteleostomi_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Euteleostomi_7') %>%
        left_join(dplyr::select(Theria_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Theria_7') %>%
        left_join(dplyr::select(Eumetazoa_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Eumetazoa_7') %>%
        left_join(dplyr::select(Homininae_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Homininae_7') %>%
        left_join(dplyr::select(Amniota_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Amniota_7') %>%
        left_join(dplyr::select(Chordata_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Chordata_sum_7') %>%
        left_join(dplyr::select(Bilateria_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Bilateria_7') %>%
        left_join(dplyr::select(Euarchontoglires_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Euarchontoglires_7') %>%
        left_join(dplyr::select(Opisthokonta_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Opisthokonta_7') %>%
        left_join(dplyr::select(Catarrhini_sum, chr.x, Average_7), by = 'chr.x') %>%
        rename('Average_7' = 'Catarrhini_7')


## Create heatmaps (chrXSpecies) for each timepoint // better visualization

rownames(Phylo_sum_0) <- Phylo_sum_0$chr.x
Phylo_sum_0 <- Phylo_sum_0[,-1]
Phylo_sum_0 <- Phylo_sum_0[-4,]
Phylo_sum_0 <- Phylo_sum_0[-7,]
Phylo_sum_0 <- Phylo_sum_0[-21,]

rownames(Phylo_sum_1) <- Phylo_sum_1$chr.x
Phylo_sum_1 <- Phylo_sum_1[,-1]
Phylo_sum_1 <- Phylo_sum_1[-4,]
Phylo_sum_1 <- Phylo_sum_1[-7,]
Phylo_sum_1 <- Phylo_sum_1[-21,]

rownames(Phylo_sum_7) <- Phylo_sum_7$chr.x
Phylo_sum_7 <- Phylo_sum_7[,-1]
Phylo_sum_7 <- Phylo_sum_7[-4,]
Phylo_sum_7 <- Phylo_sum_7[-7,]
Phylo_sum_7 <- Phylo_sum_7[-21,]

# TEST 1
# 
# H_0 <- heatmap(as.matrix(Phylo_sum_0), scale = 'column', col = heat.colors(256), 
#         main = "Species", Rowv = NA, Colv = NA)
# 
# H_1 <- heatmap(as.matrix(Phylo_sum_1), scale = 'column', col = heat.colors(256), 
#         main = "Species", Rowv = NA, Colv = NA)
# 
# H_7 <- heatmap(as.matrix(Phylo_sum_7), scale = 'column', col = heat.colors(256), 
#         main = "Species", Rowv = NA, Colv = NA)
# 
# 
# # TEST 2 - no clustering
# 
# pheatmap(as.matrix(Phylo_sum_0), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
#          fontsize_row = 5, annotation_names_col = FALSE,
#          display_numbers = TRUE, number_format = "%.2f", height = 12, width = 6)
# 
# pheatmap(as.matrix(Phylo_sum_1), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
#          fontsize_row = 5, annotation_names_col = FALSE,
#          display_numbers = TRUE, number_format = "%.2f", height = 12, width = 6)
# 
# pheatmap(as.matrix(Phylo_sum_7), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
#          fontsize_row = 5, annotation_names_col = FALSE,
#          display_numbers = TRUE, number_format = "%.2f", height = 12, width = 6)


# TEST 3 - with clustering in both rows & columns

my_colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)



col_names <- c("LUCA", "Mammalia", "H. sapiens", "Eutheria", "Eukaryota",
               "Metazoa", "Euteleostomi", "Theria", "Eumetazoa", "Homininae",
               "Amniota", "Chordata", "Bilateria", "Euarchontoglires",
               "Opisthokonta", "Catarrhini")

Phylo_sum <- data.frame(chr.x = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                                  'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                  'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                  'chr21', 'chr22', 'chrX'))



species_ordered <- c('LUCA','Eukaryota','Opisthokonta','Metazoa',
            'Eumetazoa', 'Bilateria', 'Chordata',
            'Euteleostomi', 'Amniota', 'Mammalia',
            'Theria', 'Eutheria', 'Euarchontoglires',
            'Catarrhini', 'Homininae', 'Homo Sapiens')


tiff(filename = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Heatmap/high_res_final_plot_pheatmap_02.tiff",
     width = 20, height = 15, units = "in", res = 600,
     family = "serif", compression = "lzw", bg = "white")

pheatmap(as.matrix(Phylo_sum_0),
         color=my_colors, border_color="grey90", scale="none",
         cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_method='ward.D2',
         treeheight_row = 120, treeheight_col = 120,
         cutree_rows=2,
         fontsize_row = 15, fontsize = 20, fontsize_col = 16, angle_col= "45",
         main= "0 dpi",
         annotation_names_col = FALSE,
         display_numbers = FALSE, number_format = "%.2f",
         fontsize_number = 5, number_color = 'black',
         cellheight = 25, cellwidth = 50,
         labels_col = species_ordered,
         heatmap_legend_param = list(
                 legend_height = unit(16, "cm"),  
                 grid_width = unit(1.5, "cm"),     
                 title = "", 
                 labels_gp = gpar(fontsize = 16))
)

dev.off()


tiff(filename = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Heatmap/high_res_final_plot_pheatmap_12.tiff",
     width = 20, height = 15, units = "in", res = 600,
     family = "serif", compression = "lzw", bg = "white")

pheatmap(as.matrix(Phylo_sum_1),
         color=my_colors, border_color="grey90", scale="none",
         cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_method='ward.D2',
         treeheight_row = 120, treeheight_col = 120,
         cutree_rows=2,
         fontsize_row = 15, fontsize = 20, fontsize_col = 16, angle_col= "45",
         main= "1 dpi",
         annotation_names_col = FALSE,
         display_numbers = FALSE, number_format = "%.2f",
         fontsize_number = 5, number_color = 'black',
         cellheight = 25, cellwidth = 50,
         labels_col = species_ordered,
         heatmap_legend_param = list(
                 legend_height = unit(16, "cm"),  
                 grid_width = unit(1.5, "cm"),     
                 title = "", 
                 labels_gp = gpar(fontsize = 16))
)

dev.off()


tiff(filename = "D:/Dropbox/DUTh/Thesis/SF/hgh_res/Heatmap/high_res_final_plot_pheatmap_72.tiff",
     width = 20, height = 15, units = "in", res = 600,
     family = "serif", compression = "lzw", bg = "white")

pheatmap(as.matrix(Phylo_sum_7),
         color=my_colors, border_color="grey90", scale="none",
         cluster_rows = TRUE, cluster_cols = FALSE,
         clustering_method='ward.D2',
         treeheight_row = 120, treeheight_col = 120,
         cutree_rows=2,
         fontsize_row = 15, fontsize = 20, fontsize_col = 16, angle_col= "45",
         main= "7 dpi",
         annotation_names_col = FALSE,
         display_numbers = FALSE, number_format = "%.2f",
         fontsize_number = 5, number_color = 'black',
         cellheight = 25, cellwidth = 50,
         labels_col = species_ordered,
         heatmap_legend_param = list(
                 legend_height = unit(16, "cm"),  
                 grid_width = unit(1.5, "cm"),     
                 title = "", 
                 labels_gp = gpar(fontsize = 16))
)

dev.off()
