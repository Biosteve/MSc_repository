setwd("D:/Dropbox/DUTh/Thesis/SF")

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

# Exclusive 14 genes among 9 groups

list_all <- list(
        "Set_0_low" = unlist(upset_0_low),
        "Set_1_low" = unlist(upset_1_low),
        "Set_7_low" = unlist(upset_7_low),
        "Set_0_mid" = unlist(upset_0_mid),
        "Set_1_mid" = unlist(upset_1_mid),
        "Set_7_mid" = unlist(upset_7_mid),
        "Set_0_high" = unlist(upset_0_high),
        "Set_1_high" = unlist(upset_1_high),
        "Set_7_high" = unlist(upset_7_high)
)

Upset_all <- fromList(list_all)

upset(Upset_all,
      sets = c("Set_7_high", "Set_7_mid", "Set_7_low", "Set_1_high", "Set_1_mid",
               "Set_1_low", "Set_0_high", "Set_0_mid", "Set_0_low"),
      nsets = 9,
      nintersects = NA,
      keep.order = TRUE,
      show.numbers = 'Yes',
      scale.intersections = 'identity', 
      scale.sets = 'identity')

all_genes <- unique(unlist(list_all))

all_genes <-as.data.frame(all_genes)
all_genes <- all_genes %>% mutate(ID_number = row_number()) %>%
        relocate(ID_number, all_genes)
all_genes$ID_number <- as.character(all_genes$ID_number)

gene_name_list <- data.frame(gene_name = all_genes, stringsAsFactors = FALSE)

#==============================================================================

exclusive14_among_9_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
        Upset_all$Set_0_mid == 1 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 1 &
        Upset_all$Set_1_mid == 1 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 1 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 1, 
]

gene_ID_14 <- row.names(exclusive14_among_9_groups)
gene_14 <- data.frame(gene_ID_14)
colnames(gene_14)[colnames(gene_14) == "gene_ID_14"] <- 'ID_number'
gene_14 <- gene_14 %>% left_join(all_genes, by='ID_number')

colnames(gene_14)[colnames(gene_14) == 'all_genes'] <- 'Gene'
gene_14_list <- merge(x=gene_14, y=Average_0, by ='Gene' )
gene_14_list_sum <- gene_14_list %>% count(Gene)

##RESULT : based on two different methods, the 14 genes are the same!!!

#=============================================================================

# Exclusive 6 genes among 8 groups

exclusive6_among_8_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
        Upset_all$Set_0_mid == 1 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 1 &
        Upset_all$Set_1_mid == 1 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 0 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 1, 
]

gene_ID_6 <- row.names(exclusive6_among_8_groups)
gene_6 <- data.frame(gene_ID_6)
colnames(gene_6)[colnames(gene_6) == "gene_ID_6"] <- 'ID_number'
gene_6 <- gene_6 %>% left_join(all_genes, by='ID_number')

colnames(gene_6)[colnames(gene_6) == 'all_genes'] <- 'Gene'
gene_6_list <- merge(x=gene_6, y=Average_0, by = 'Gene')
gene_6_list_sum <- gene_6_list %>% count(Gene)

#=============================================================================

# Exclusive 7 genes among 6 groups

exclusive7_among_6_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 0 &
        Upset_all$Set_1_mid == 1 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 1 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_7a <- row.names(exclusive7_among_6_groups)
gene_7a <- data.frame(gene_ID_7a)
colnames(gene_7a)[colnames(gene_7a) == "gene_ID_7a"] <- 'ID_number'
gene_7a <- gene_7a %>% left_join(all_genes, by='ID_number')

colnames(gene_7a)[colnames(gene_7a) == 'all_genes'] <- 'Gene'
gene_7a_list <- merge(x=gene_7a, y=Average_0, by = 'Gene')
gene_7a_list_sum <- gene_7a_list %>% count(Gene)

#=============================================================================

# Exclusive 7 genes among 4 groups

exclusive7_among_4_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 0 &
        Upset_all$Set_1_mid == 1 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 0 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_7b <- row.names(exclusive7_among_4_groups)
gene_7b <- data.frame(gene_ID_7b)
colnames(gene_7b)[colnames(gene_7b) == "gene_ID_7b"] <- 'ID_number'
gene_7b <- gene_7b %>% left_join(all_genes, by='ID_number')

colnames(gene_7b)[colnames(gene_7b) == 'all_genes'] <- 'Gene'
gene_7b_list <- merge(x=gene_7b, y=Average_0, by = 'Gene')
gene_7b_list_sum <- gene_7b_list %>% count(Gene)

#=============================================================================

# Exclusive 799 genes among 3 groups

exclusive799_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 0 &
        Upset_all$Set_1_low == 1 &
        Upset_all$Set_1_mid == 0 &
        Upset_all$Set_1_high == 0 &
        Upset_all$Set_7_low == 1 &
        Upset_all$Set_7_mid == 0 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_799 <- row.names(exclusive799_among_3_groups)
gene_799 <- data.frame(gene_ID_799)
colnames(gene_799)[colnames(gene_799) == "gene_ID_799"] <- 'ID_number'
gene_799 <- gene_799 %>% left_join(all_genes, by='ID_number')

colnames(gene_799)[colnames(gene_799) == 'all_genes'] <- 'Gene'
gene_799_list <- merge(x=gene_799, y=Average_0, by = 'Gene')
gene_799_list_sum <- gene_799_list %>% count(Gene)

#=============================================================================

# Exclusive 1305 genes among 3 groups

exclusive1305_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 0 &
        Upset_all$Set_1_low == 0 &
        Upset_all$Set_1_mid == 0 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 1 &
        Upset_all$Set_7_mid == 0 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_1305 <- row.names(exclusive1305_among_3_groups)
gene_1305 <- data.frame(gene_ID_1305)
colnames(gene_1305)[colnames(gene_1305) == "gene_ID_1305"] <- 'ID_number'
gene_1305 <- gene_1305 %>% left_join(all_genes, by='ID_number')

colnames(gene_1305)[colnames(gene_1305) == 'all_genes'] <- 'Gene'
gene_1305_list <- merge(x=gene_1305, y=Average_0, by = 'Gene')
gene_1305_list_sum <- gene_1305_list %>% count(Gene)

#=============================================================================

# Exclusive 10337 genes among 3 groups

exclusive10337_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 0 &
        Upset_all$Set_1_mid == 1 &
        Upset_all$Set_1_high == 0 &
        Upset_all$Set_7_low == 0 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_10337 <- row.names(exclusive10337_among_3_groups)
gene_10337 <- data.frame(gene_ID_10337)
colnames(gene_10337)[colnames(gene_10337) == "gene_ID_10337"] <- 'ID_number'
gene_10337 <- gene_10337 %>% left_join(all_genes, by='ID_number')

colnames(gene_10337)[colnames(gene_10337) == 'all_genes'] <- 'Gene'
gene_10337_list <- merge(x=gene_10337, y=Average_0, by = 'Gene')
gene_10337_list_sum <- gene_10337_list %>% count(Gene)

#=============================================================================

# Exclusive 13771 genes among 3 groups

exclusive13771_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
        Upset_all$Set_0_mid == 0 &
        Upset_all$Set_0_high == 1 &
        Upset_all$Set_1_low == 0 &
        Upset_all$Set_1_mid == 0 &
        Upset_all$Set_1_high == 1 &
        Upset_all$Set_7_low == 0 &
        Upset_all$Set_7_mid == 1 &
        Upset_all$Set_7_high == 0, 
]

gene_ID_13771 <- row.names(exclusive13771_among_3_groups)
gene_13771 <- data.frame(gene_ID_13771)
colnames(gene_13771)[colnames(gene_13771) == "gene_ID_13771"] <- 'ID_number'
gene_13771 <- gene_13771 %>% left_join(all_genes, by='ID_number')

colnames(gene_13771)[colnames(gene_13771) == 'all_genes'] <- 'Gene'
gene_13771_list <- merge(x=gene_13771, y=Average_0, by = 'Gene')
gene_13771_list_sum <- gene_13771_list %>% count(Gene)

#=============================================================================

# Exclusive 5266 genes among 3 groups

exclusive5266_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
                Upset_all$Set_0_mid == 1 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 0 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 1 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 1, 
]

gene_ID_5266 <- row.names(exclusive5266_among_3_groups)
gene_5266 <- data.frame(gene_ID_5266)
colnames(gene_5266)[colnames(gene_5266) == "gene_ID_5266"] <- 'ID_number'
gene_5266 <- gene_5266 %>% left_join(all_genes, by='ID_number')

colnames(gene_5266)[colnames(gene_5266) == 'all_genes'] <- 'Gene'
gene_5266_list <- merge(x=gene_5266, y=Average_0, by = 'Gene')
gene_5266_list_sum <- gene_5266_list %>% count(Gene)

#=============================================================================

# Exclusive 4029 genes among 3 groups

exclusive4029_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
                Upset_all$Set_0_mid == 0 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 1 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 0 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 1, 
]

gene_ID_4029 <- row.names(exclusive4029_among_3_groups)
gene_4029 <- data.frame(gene_ID_4029)
colnames(gene_4029)[colnames(gene_4029) == "gene_ID_4029"] <- 'ID_number'
gene_4029 <- gene_4029 %>% left_join(all_genes, by='ID_number')

colnames(gene_4029)[colnames(gene_4029) == 'all_genes'] <- 'Gene'
gene_4029_list <- merge(x=gene_4029, y=Average_0, by = 'Gene')
gene_4029_list_sum <- gene_4029_list %>% count(Gene)

#=============================================================================

# Exclusive 3219 genes among 3 groups

exclusive3219_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
                Upset_all$Set_0_mid == 0 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 0 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 1 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 1, 
]

gene_ID_3219 <- row.names(exclusive3219_among_3_groups)
gene_3219 <- data.frame(gene_ID_3219)
colnames(gene_3219)[colnames(gene_3219) == "gene_ID_3219"] <- 'ID_number'
gene_3219 <- gene_3219 %>% left_join(all_genes, by='ID_number')

colnames(gene_3219)[colnames(gene_3219) == 'all_genes'] <- 'Gene'
gene_3219_list <- merge(x=gene_3219, y=Average_0, by = 'Gene')
gene_3219_list_sum <- gene_3219_list %>% count(Gene)

#=============================================================================

# Exclusive 2931 genes among 3 groups

exclusive2931_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
                Upset_all$Set_0_mid == 0 &
                Upset_all$Set_0_high == 1 &
                Upset_all$Set_1_low == 0 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 1 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 1, 
]

gene_ID_2931 <- row.names(exclusive2931_among_3_groups)
gene_2931 <- data.frame(gene_ID_2931)
colnames(gene_2931)[colnames(gene_2931) == "gene_ID_2931"] <- 'ID_number'
gene_2931 <- gene_2931 %>% left_join(all_genes, by='ID_number')

colnames(gene_2931)[colnames(gene_2931) == 'all_genes'] <- 'Gene'
gene_2931_list <- merge(x=gene_2931, y=Average_0, by = 'Gene')
gene_2931_list_sum <- gene_2931_list %>% count(Gene)

#=============================================================================

# Exclusive 2169a genes among 3 groups

exclusive2169a_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
                Upset_all$Set_0_mid == 1 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 1 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 0 &
                Upset_all$Set_7_low == 1 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 0, 
]

gene_ID_2169a <- row.names(exclusive2169a_among_3_groups)
gene_2169a <- data.frame(gene_ID_2169a)
colnames(gene_2169a)[colnames(gene_2169a) == "gene_ID_2169a"] <- 'ID_number'
gene_2169a <- gene_2169a %>% left_join(all_genes, by='ID_number')

colnames(gene_2169a)[colnames(gene_2169a) == 'all_genes'] <- 'Gene'
gene_2169a_list <- merge(x=gene_2169a, y=Average_0, by = 'Gene')
gene_2169a_list_sum <- gene_2169a_list %>% count(Gene)

#=============================================================================

# Exclusive 2169b genes among 3 groups

exclusive2169b_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
                Upset_all$Set_0_mid == 1 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 1 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 0 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 0 &
                Upset_all$Set_7_high == 1, 
]

gene_ID_2169b <- row.names(exclusive2169b_among_3_groups)
gene_2169b <- data.frame(gene_ID_2169b)
colnames(gene_2169b)[colnames(gene_2169b) == "gene_ID_2169b"] <- 'ID_number'
gene_2169b <- gene_2169b %>% left_join(all_genes, by='ID_number')

colnames(gene_2169b)[colnames(gene_2169b) == 'all_genes'] <- 'Gene'
gene_2169b_list <- merge(x=gene_2169b, y=Average_0, by = 'Gene')
gene_2169b_list_sum <- gene_2169b_list %>% count(Gene)

#=============================================================================

# Exclusive 2106 genes among 3 groups

exclusive2106_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 0 & 
                Upset_all$Set_0_mid == 1 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 0 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 1 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 1 &
                Upset_all$Set_7_high == 0, 
]

gene_ID_2106 <- row.names(exclusive2106_among_3_groups)
gene_2106 <- data.frame(gene_ID_2106)
colnames(gene_2106)[colnames(gene_2106) == "gene_ID_2106"] <- 'ID_number'
gene_2106 <- gene_2106 %>% left_join(all_genes, by='ID_number')

colnames(gene_2106)[colnames(gene_2106) == 'all_genes'] <- 'Gene'
gene_2106_list <- merge(x=gene_2106, y=Average_0, by = 'Gene')
gene_2106_list_sum <- gene_2106_list %>% count(Gene)

#=============================================================================

# Exclusive 1354 genes among 3 groups

exclusive1354_among_3_groups <- Upset_all[
        Upset_all$Set_0_low == 1 & 
                Upset_all$Set_0_mid == 0 &
                Upset_all$Set_0_high == 0 &
                Upset_all$Set_1_low == 0 &
                Upset_all$Set_1_mid == 0 &
                Upset_all$Set_1_high == 1 &
                Upset_all$Set_7_low == 0 &
                Upset_all$Set_7_mid == 1 &
                Upset_all$Set_7_high == 0, 
]

gene_ID_1354 <- row.names(exclusive1354_among_3_groups)
gene_1354 <- data.frame(gene_ID_1354)
colnames(gene_1354)[colnames(gene_1354) == "gene_ID_1354"] <- 'ID_number'
gene_1354 <- gene_1354 %>% left_join(all_genes, by='ID_number')

colnames(gene_1354)[colnames(gene_1354) == 'all_genes'] <- 'Gene'
gene_1354_list <- merge(x=gene_1354, y=Average_0, by = 'Gene')
gene_1354_list_sum <- gene_1354_list %>% count(Gene)

#=============================================================================


