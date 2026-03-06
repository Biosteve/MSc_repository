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


##Data_0##
data_0<-read.table("d0_model_ladatlas_vs_allGenes.bed")

summary(data_0)
test_0<-data_0[1:8]
sum(is.na(test_0)) 

# test_0<-test_0 %>% mutate(newV6=V6**2)
# test_0<-test_0 %>% mutate(newV7=V7**2)
# test_0<-test_0 %>% mutate(newV8=V8**2)
# test_0<-test_0 %>% mutate(sum=newV6+newV7+newV8)
# test_0<-test_0 %>% mutate(dist=sqrt(sum))

# to make it in a for loop! // να υπολογίζει αυτόματα για V6-V8 και μετά να υπολογίζει το sqrt


calc_Dist <- function(test_0){
        test_0<-test_0 %>% mutate(sum=V6**2+V7**2+V8**2)
        test_0<-test_0 %>% mutate(dist=sqrt(sum))
        test_0<-test_0 %>% mutate(size=V3-V2)
        colnames(test_0)=c("chr", "startPOS", "endPOS", "Gene", "Sis", "X", "Y", "Z", "sum", "Dist", "Gene size")
        return(test_0)
}

Dist_0=calc_Dist(test_0)
head(Dist_0)
summary(Dist_0)


Subset_sis<-function(x){
        test_A<-x[which(x$"Sis"=="A"), ]
        test_B<-x[which(x$"Sis"=="B"), ]
        return(list("A"=test_A, "B"=test_B))
}

SisA_0=Subset_sis(Dist_0)$A
SisB_0=Subset_sis(Dist_0)$B
summary(SisA_0$Dist)
summary(SisB_0$Dist)



##Data_1##

data_1<-read.table("d1_model_ladatlas_vs_allGenes.bed")
summary(data_1)
test_1<-data_1[1:8]
summary(test_1)

# Checking for missing values, remove them and cross-check (from 116242 to 113856, removing 2386 rows)
sum(is.na(test_1))
test_1_noNA<-na.omit(test_1)
sum(is.na(test_1_noNA))
summary(test_1_noNA)



Dist_1=calc_Dist(test_1)
head(Dist_1)
summary(Dist_1)


SisA_1=Subset_sis(Dist_1)$A
SisB_1=Subset_sis(Dist_1)$B
summary(SisA_1$Dist)
summary(SisB_1$Dist)


##Data_7##

data_7<-read.table("d7_model_ladatlas_vs_allGenes.bed")
summary(data_7)
test_7<-data_7[1:8]
summary(test_7)


# Checking for missing values, remove them and cross-check (from 116242 to 111357, removing 4885 rows)
sum(is.na(test_7))
test_7_noNA<-na.omit(test_7)
sum(is.na(test_7_noNA))
summary(test_7_noNA)



Dist_7=calc_Dist(test_7)
head(Dist_7)
summary(Dist_7)
SisA_7=Subset_sis(Dist_7)$A
SisB_7=Subset_sis(Dist_7)$B

summary(SisA_7$Dist)
summary(SisB_7$Dist)
head(SisB_7)


##GGpolt Graphs##

graphing_0<- ggplot(Dist_0, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="0 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

graphing_1<- ggplot(Dist_1, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="1 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

graphing_7<- ggplot(Dist_7, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,6)) + 
        geom_boxplot(alpha=1) + labs(title="7 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

ggarrange(graphing_0, graphing_1, graphing_7, ncol =3 , nrow = 1, common.legend = TRUE, legend='bottom')

##Create datasets that contain only NAs and save them (need the function to do it)

OnlyNA_d1<- test_1 %>% filter_all(any_vars(is.na(.)))
OnlyNA_d7<- test_7 %>% filter_all(any_vars(is.na(.)))


##Using unique, I can examine what values are contained within this datasets

unique(OnlyNA_d1$V1) #--> chr8
unique(OnlyNA_d1$V5) #--> A
unique(OnlyNA_d7$V1) #--> chr4, chrX
unique(OnlyNA_d7$V5) #--> B

##What percentage of the total values are the NAs?

D1_chr8<-subset(test_1, V1=='chr8') ## creates a dataset with only the chr8
table(D1_chr8$V5) # 2386 rows (50%) are missing values and only SiS_A
D7_chr4<-subset(test_7, V1=='chr4') ## creates a dataset with only the chr4
table(D7_chr4$V5) # 2510 rows (50%) are missing values and only SiS_B
D7_chrX<-subset(test_7, V1=='chrX') ## creates a dataset with only the chrX
table(D7_chrX$V5) # 2375 rows (50%) are missing values and only SiS_B

##How to merge common data from the 3 data_sets?? combine <-cbind(d0, d1, d7)
 #If exclude all the Chr data that have NAs, then we will work with 101,700 genes
 # and we will go w/o chr4,8,X!!!
 # function that creates a new dataset w/o the chr4,8,X!!

# Dist_0_NA<-subset(Dist_0, chr!= 'chrX')
# Dist_0_NA<-subset(Dist_0_NA, chr!= 'chr4')
# Dist_0_NA<-subset(Dist_0_NA, chr!= 'chr8')
# 
# unique(Dist_0_NA$chr) ## You need 101700 rows // Check whether the chr are excluded
# 
# Dist_1_NA<-subset(Dist_1, chr!= 'chrX')
# Dist_1_NA<-subset(Dist_1_NA, chr!= 'chr4')
# Dist_1_NA<-subset(Dist_1_NA, chr!= 'chr8')
# 
# unique(Dist_1_NA$chr) ## You need 101700 rows // Check whether the chr are excluded
# 
# Dist_7_NA<-subset(Dist_7, chr!= 'chrX')
# Dist_7_NA<-subset(Dist_7_NA, chr!= 'chr4')
# Dist_7_NA<-subset(Dist_7_NA, chr!= 'chr8')


unique(Dist_7_NA$chr) ## You need 101700 rows // Check whether the chr are excluded

##create a function that creates the Dist_X_NA datasets!

Dist_NA <- function(Dist_0){
        Dist_a<-subset(Dist_0, chr!='chrX')
        Dist_b<-subset(Dist_a, chr!='chr4')
        Dist_0<-subset(Dist_b, chr!='chr8')
        return(Dist_0)
}


Dist_0_NA=Dist_NA(Dist_0)
head(Dist_0_NA)
summary(Dist_0_NA)

Dist_1_NA=Dist_NA(Dist_1)
head(Dist_1_NA)
summary(Dist_1_NA)

Dist_7_NA=Dist_NA(Dist_7)
head(Dist_7_NA)
summary(Dist_7_NA)

SisA_0_NA=Subset_sis(Dist_0_NA)$A
SisB_0_NA=Subset_sis(Dist_0_NA)$B
summary(SisA_0_NA$Dist)
summary(SisB_0_NA$Dist)

SisA_1_NA=Subset_sis(Dist_1_NA)$A
SisB_1_NA=Subset_sis(Dist_1_NA)$B
summary(SisA_1_NA$Dist)
summary(SisB_1_NA$Dist)

SisA_7_NA=Subset_sis(Dist_7_NA)$A
SisB_7_NA=Subset_sis(Dist_7_NA)$B
summary(SisA_7_NA$Dist)
summary(SisB_7_NA$Dist)

# graphs with only the common chrs in all 3 datasets

graphing_0_NA<- ggplot(Dist_0_NA, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) + 
        stat_compare_means(comparisons = list(c('A', 'B')),
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="0 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5),
              axis.text.x = element_text(size = 15, color = 'black'),
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              axis.text.y = element_text(size = 15, color = 'black'),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

graphing_1_NA<- ggplot(Dist_1_NA, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) + 
        stat_compare_means(comparisons = list(c('A', 'B')),
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="1 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5),
              axis.text.x = element_text(size = 15, color = 'black'),
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              axis.text.y = element_text(size = 15, color = 'black'),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

graphing_7_NA<- ggplot(Dist_7_NA, aes(x=Sis, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) + 
        stat_compare_means(comparisons = list(c('A', 'B')),
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="7 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_text(size = 20, vjust=-1, hjust=0.5),
              axis.text.x = element_text(size = 15, color = 'black'),
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_text(size = 20, vjust=2, hjust=0.5),
              axis.text.y = element_text(size = 15, color = 'black'),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

final_graph_noNA <- ggarrange(graphing_0_NA, graphing_1_NA, graphing_7_NA,
                              ncol =3 , nrow = 1, common.legend = FALSE)

ggsave("high_res_final_graph_noNA.tiff",
       plot = final_graph_noNA,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# graphs based on chr

correct_order <- c("chr1", "chr2", "chr3", "chr5", "chr6", "chr7", "chr9", "chr10", 
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
                   "chr18", "chr19", "chr20", "chr21", "chr22")

Dist_0_NA$chr <- factor(Dist_0_NA$chr, levels = correct_order)

graphing_0_NA_chr<- ggplot(Dist_0_NA, aes(x=chr, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) +  
        stat_compare_means(aes(group=Sis),
                           method='wilcox.test',
                           label = 'p.signif') +
        labs(title="0 dpi", x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(),
              axis.text.x= element_text(family='serif'),
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

Dist_1_NA$chr <- factor(Dist_1_NA$chr, levels = correct_order)

graphing_1_NA_chr<- ggplot(Dist_1_NA, aes(x=chr, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) + 
        stat_compare_means(aes(group=Sis),
                           method='wilcox.test',
                           label = 'p.signif',
                           label.y=6) +
        labs(title="1 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())

Dist_7_NA$chr <- factor(Dist_7_NA$chr, levels = correct_order)

graphing_7_NA_chr<- ggplot(Dist_7_NA, aes(x=chr, y=Dist, fill=Sis)) + 
        scale_y_continuous(limits = c(0,7)) + 
        geom_boxplot(alpha=1) + 
        stat_compare_means(aes(group=Sis),
                           method='wilcox.test',
                           label = 'p.signif',
                           label.y=6) +
        labs(title="7 dpi") + 
        labs(x= "Sister chromatids", y="Distance from nuclear center") + 
        theme(axis.line.x = element_line(linewidth = 1, colour = "black"), 
              axis.title.x = element_blank(), 
              axis.line.y = element_line(linewidth = 1, colour = "black"), 
              axis.title.y = element_blank(),
              legend.position="bottom", legend.direction="horizontal", 
              legend.title = element_blank(), 
              plot.title = element_text(family = "serif", hjust=0.5), 
              text = element_text(family = "serif", face='bold'), 
              panel.background = element_blank())


# graphing_0_NA_chr<- ggplot(Dist_0_NA, aes(x=chr, y=Dist, color=Sis)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Mean_Sis_Dist_0_NA") + labs(x= "Sis", y="Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_1_NA_chr<- ggplot(Dist_1_NA, aes(x=chr, y=Dist, color=Sis)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Mean_Sis_Dist_1_NA") + labs(x= "Sis", y="Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))
# graphing_7_NA_chr<- ggplot(Dist_7_NA, aes(x=chr, y=Dist, color=Sis)) + scale_y_continuous(limits = c(0,6)) + geom_boxplot(alpha=1) + labs(title="Mean_Sis_Dist_7_NA") + labs(x= "Sis", y="Distance from Center") + theme(axis.line.x = element_line(linewidth = 0.5, colour = "black"), legend.position="bottom", legend.direction="horizontal", legend.title = element_blank(), plot.title = element_text(family = "Calibri"), text = element_text(family = "Calibri"))

combined_plot_Dist_NA <- ggarrange(graphing_0_NA_chr, graphing_1_NA_chr, graphing_7_NA_chr,ncol = 1, nrow = 3, common.legend = TRUE, legend='bottom')

final_plot_Dist_NA <- annotate_figure(
        combined_plot_Dist_NA,
        bottom = textGrob('Chromosomes', gp=gpar(fontfamily='serif', fontface='bold', fontsize=20)),
        left = textGrob("Distances from nuclear center", rot=90, gp=gpar(fontfamily='serif', fontface='bold', fontsize=20))
)

print(final_plot_Dist_NA)

ggsave("high_res_final_plot_Dist_NA.tiff",
       plot = final_plot_Dist_NA,
       path = "D:/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8,units = 'in', dpi = 900, device = 'tiff',
       compression = 'lzw', bg='white')


# To examine whether genes have different behavior? Check genes based on chr?

# ##Dataset_0 #Does the data follow normal distribution.??
# 
# Test_0_SisA<- which(Dist_0_NA$Sis=="A")
# Test_0_SisB<- which(Dist_0_NA$Sis=="B")
# 
# t.test(Dist_0_NA$Dist[Test_0_SisA], Dist_0_NA$Dist[Test_0_SisB], conf.level = 0.95)
# 
# ##Dataset_1
# 
# Test_1_SisA<- which(Dist_1_NA$Sis=="A")
# Test_1_SisB<- which(Dist_1_NA$Sis=="B")
# 
# t.test(Dist_1_NA$Dist[Test_1_SisA], Dist_1_NA$Dist[Test_1_SisB], conf.level = 0.95)
# 
# ##Dataset_7
# 
# Test_7_SisA<- which(Dist_7_NA$Sis=="A")
# Test_7_SisB<- which(Dist_7_NA$Sis=="B")
# 
# t.test(Dist_7_NA$Dist[Test_7_SisA], Dist_7_NA$Dist[Test_7_SisB], conf.level = 0.95)

# Human Genome Assembly GRCh38.p14 - For Chromosome lengths

# Creating a new dataset to calculate the average value of each gene
# between the two sis for each timepoint

# AverageSIS.0 <- cbind(SisA_0, SisB_0)
# AverageSIS_A0 <- AverageSIS.0[-c(6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 22)]
# AverageSIS_0 <- AverageSIS_A0 %>% mutate(mean=( Dist+Dist.1)/2)
# AverageSIS_0 <- AverageSIS_0 %>% rename(Sis.A = Sis, Sum.A = sum, Dist.A = Dist, Sis.B = Sis.1, Sum.B = sum.1, Dist.B = Dist.1)
# 
# 
# AverageSIS.1 <- cbind(SisA_1, SisB_1)
# AverageSIS_A1 <- AverageSIS.1[-c(6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 22)]
# AverageSIS_1 <- AverageSIS_A1 %>% mutate(mean=( Dist+Dist.1)/2)
# AverageSIS_1 <- AverageSIS_1 %>% rename(Sis.A = Sis, Sum.A = sum, Dist.A = Dist, Sis.B = Sis.1, Sum.B = sum.1, Dist.B = Dist.1)
# 
# 
# AverageSIS.7 <- cbind(SisA_7, SisB_7)
# AverageSIS_A7 <- AverageSIS.7[-c(6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 22)]
# AverageSIS_7 <- AverageSIS_A7 %>% mutate(mean=( Dist+Dist.1)/2)
# AverageSIS_7 <- AverageSIS_7 %>% rename(Sis.A = Sis, Sum.A = sum, Dist.A = Dist, Sis.B = Sis.1, Sum.B = sum.1, Dist.B = Dist.1)



# To create a function that does all the above

Average <-function(x, y, z, w){
        x <- cbind(y, z)
        x.1 <- x[-c(6, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 22)]
        colnames(x.1)=c("chr", "startPOS", "endPOS", "Gene", "Sis.A", "sum.A", "Dist.A", "Sis.B", "sum.B", "Dist.B")
        x.2 <- x.1 %>% mutate(average=(Dist.A+Dist.B)/2)
        x.3 <- x.2 %>% add_column(timepoint= w, .before = "average")
        return(x.3)
        
}

Average_0 = Average(Dist_0_NA, SisA_0_NA, SisB_0_NA, "d0")
Average_1 = Average(Dist_1_NA, SisA_1_NA, SisB_1_NA, "d1")
Average_7 = Average(Dist_7_NA, SisA_7_NA, SisB_7_NA, "d7") 


# remove the data that appear more than 4 times in Gene column in Average_0,1,7
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

# Then reassign the Average groups with the single copies genes

Average_0 <- filtered_Average_0
Average_1 <- filtered_Average_1
Average_7 <- filtered_Average_7

# Average_All <- rbind(Average_0, Average_1, Average_7)

#Pivot.of.average <- Average_All %>% group_by(chr) %>% summarise(Average.0 = mean(average), SD.0 = sd(average),Average.1 = mean(average.1), SD.1 = sd(average.1), Average.7 = mean(average.7), SD.7 = sd(average.7))

# ggplot(data=Average_All) + geom_boxplot(aes(timepoint, average, color=timepoint)) +
#         scale_fill_manual(values = T) +
#         facet_wrap(~chr) +
#         theme_minimal()
# 
# Average_sub <- subset(Average_All, select = c(chr, Dist.A, Dist.B, timepoint))
# Average_sub <- na.omit(Average_sub)

Sum_d0 <- summarise(Average_0, .by = chr, mean.SiSa = mean(Dist.A), mean.SiSb = mean(Dist.B), timepoint)
Sum_d1 <- summarise(Average_1, .by = chr, mean.SiSa = mean(Dist.A), mean.SiSb = mean(Dist.B), timepoint)
Sum_d7 <- summarise(Average_7, .by = chr, mean.SiSa = mean(Dist.A), mean.SiSb = mean(Dist.B), timepoint)

Sum_All <- cbind (Sum_d0, Sum_d1, Sum_d7)
colnames(Sum_All)=c("chr", "SisA.0", "SisB.0", "d.0", "chr.1", "SisA.1", "SisB.1", "d.1", "chr.7", "SisA.7", "SisB.7", "d.7")

Sum_All <- na.omit(Sum_All)

Sum_All <- Sum_All[, !names(Sum_All) %in% c("chr.1", "chr.7")]


##Calculate the differences between the two Sis

Group_All <- Sum_All %>% group_by(chr) %>%
        summarise(across(c(SisA.0, SisB.0, SisA.1, SisB.1, SisA.7, SisB.7), mean), .groups = 'drop') %>%
        mutate(Average.0 = (SisA.0+SisB.0)/2, Average.1 = (SisA.1+SisB.1)/2, Average.7 = (SisA.7+SisB.7)/2, Diff_0 = abs(SisA.0 - SisB.0), Diff_1 = abs(SisA.1 - SisB.1), Diff_7 = abs(SisA.7 - SisB.7)) %>%
        as.data.frame()

ALL <- merge(Group_All, Genome_density, by='chr')

ALL <- ALL[, -c(2, 3, 4, 5, 6, 7, 15, 16)]
    
write_xlsx(ALL, 'ALL_data.xlsx')
