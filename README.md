# MSc_repository (Longitudinal Systems Biology Analysis Pipeline)
The code that was created and used for the results of my MSc thesis titled "Analysis of the evolutionary history of eukaryotic genes in relation to their spatial arrangement in the eukaryotic nucleus"


📌 Project Overview
* An end-to-end data pipeline built in R to process and model longitudinal genomic data, analyzing over 350,000 observations across 3 timepoints.

💻 Tech Stack
* Language: R
* Data Wrangling: Tidyverse (dplyr, tidyr)
* Bioinformatics & Systems Biology: Bioconductor, ComplexHeatmap, GenomicRanges, UpSetR
* Visualization: ggplot2


⚙️ Methodology & Pipeline
* Data Ingestion & Cleaning: Noise reduction and data standardization on a 116k+ row dataset.
* Feature Engineering & Statistics: Feature extraction and variance calculation.
* Network Modeling: Identification of hidden regulatory patterns in a high-dimensional environment.

📊 Key Visualizations :arrow_right: [HeatMap](output/Pheatmap_0days.tiff)

🚀 How to Run
- To run the pipeline locally: Execute scripts in this order:
- 
1. Radia_positioning.R
2. Gene density.R
3. Phylostratigraphy.R
4. UpSet_plots.R
5. Upset_exclusive_filtered.R
6. gProfPlot2.R

:bangbang:Data availability:bangbang:
- The data that support the findings of this Thesis are available upon reasonable request
- Sample datasets are created for reproducible purposes 
