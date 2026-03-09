
##Data_read##

ref <- read.table("hg38.genome")
ref <- ref[-(24:455),]
colnames(ref)=c('chr', 'reference')

ref2 <- data.frame(chr=c('chr21'),
                   reference = c('46709983'))
ref <- rbind (ref, ref2)


data_read<-read.table("merged.bed")
summary(data_read)
colnames(data_read)=c("chr", "startPOS", "endPOS")


size <- data_read %>%
        mutate(size = endPOS - startPOS)

Genome_size <- size %>% group_by(chr) %>%
        summarise(chr_len = sum(size))

Genome_size <- Genome_size %>% left_join(ref, by = 'chr')

Genome_density <- Genome_size %>%
        mutate(chr_len = as.numeric(as.character(chr_len)),
               reference = as.numeric(as.character(reference)),
               gen.den = (chr_len / reference)*100)

# function to format the labels as `10^x expressions
scientific_10 <- function(x){
        coeff <- x / 1e8
        parse(text=paste0(coeff, " %*% 10^8"))
                          }

#Correlation gene density - chromosome length

plot_data <- Genome_density %>%
  filter(!chr %in% c('chr4', 'chr8', 'chrX')) %>%
  mutate(
      v_adj = case_when(
        chr == 'chr2' ~ -1.5,
        chr == 'chr12' ~ 1,
        chr == 'chr11' ~ 1,
        chr == 'chr13' ~ 1,
        chr %in% c('chr8', 'chr6', 'chr17') ~ 2,
        TRUE ~ 2.5
    ),
    h_adj = case_when(
      chr == 'chr5' ~ -0.5,
      chr == 'chr11' ~ -0.5,
      chr == 'chr12' ~ 1.5,
      chr == 'chr13' ~ 1.5,
      TRUE ~ 0.5
    )
  )

Gen_den_plot <- ggplot(data = plot_data, 
                       aes(x = reference, y = gen.den)) +
  geom_point(fill = 'darkred', color = 'black', size = 8, shape = 21, stroke=1.2) +
  geom_text(aes(label = chr,
                vjust=v_adj,
            hjust=h_adj,
                family = 'serif',
                fontface='bold',
                size=6)) +
  scale_x_continuous(labels = function(x) sprintf("%.1f", x/1e8)) +
  labs(y="Gene density (%)", x=expression(bold("Chromosome size ("*10^8* " bp)"))) +
  theme_classic() +
  theme(axis.line.x = element_line(linewidth = 1, colour = "black"),
        axis.text.x = element_text(size = 15, color = 'black'),
        axis.title.x = element_text(size = 20, vjust=1, hjust=0.5),
        axis.line.y = element_line(linewidth = 1, colour = "black"),
        axis.text.y = element_text(size = 15, color = 'black'),
        axis.title.y = element_text(size = 20, vjust=1, hjust=0.5),
        legend.position="bottom", legend.direction="horizontal",
        legend.title = element_blank(),
        plot.title = element_text(family = "serif", hjust=0.5),
        text = element_text(family = "serif", face='bold'),
        panel.background = element_blank(),)


print(Gen_den_plot)


# Saves in the directory
ggsave("high_res_gen_den_plot.tiff",
       plot = Gen_den_plot,
       path = "C:/Users/User/Dropbox/DUTh/Thesis/SF/hgh_res",
       width = 15, height = 8, 
       units = 'in', dpi = 900, device = 'tiff', compression = 'lzw')
