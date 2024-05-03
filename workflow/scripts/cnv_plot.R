log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))

# cnv <- read.delim("results/dataset/copy_number_variants_dataset.tsv", sep = "\t", header = TRUE)
# repeats <- read.delim("results/references/repeats.bed", sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"))
# metadata <- read.delim("config/sample_metadata.csv", sep = ",", header = TRUE)
# chrom_names <- read.delim("config/chromosome_names.csv", sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))
print("Reading files")
cnv <- read.delim(snakemake@input[[1]], sep = "\t", header = TRUE)
repeats <- read.delim(snakemake@input[[2]], sep = "\t", header = FALSE, col.names = c("Accession", "Start", "End", "Repeat_type"))
metadata <- read.delim(snakemake@input[[3]], sep = ",", header = TRUE)
chrom_names <- read.delim(snakemake@input[[4]], sep = ",", header = FALSE, col.names = c("Lineage", "Accession", "Chromosome"))

metadata <- metadata %>% 
    select(sample, lineage = group)

print("Joining files")
cnv <- left_join(cnv, metadata, by = c("Sample" = "sample"))
repeats <- left_join(repeats, chrom_names, by = "Accession")
repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))

print("Separating lineages")
lineage <- unique(cnv$lineage)
for (lin in lineage){
    df_lin <- cnv %>% filter(lineage == lin)
    df_repeats <- repeats %>% filter(Lineage == lin)
    assign(lin, df_lin)
    assign(paste0(lin, "_repeats"), df_repeats)
}

set2 <- rev(brewer.pal(8, "Set2")[1:6])
s_colors <- set2[1:length(unique(cnv$Structure))]

print("Plotting")
for (lin in lineage){
    print(lin)
    print("Plotting CNV")
    plot_cnv <- ggplot(get(lin)) +
        ylim(1,1)+
        geom_segment(aes(x=Start, xend=End, y = 1, yend = 1 , color = Structure), linewidth = 5)+
        scale_color_manual(values = s_colors)+
        facet_grid(Sample~Accession, scales = "free")+
        scale_x_continuous(breaks = seq(0, max(get(lin)$End, na.rm = TRUE), by = 500000), labels = scales::label_comma())+
        theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 0.5),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            strip.text.y = element_text(angle = 0),
            strip.text.x = element_blank(),
            panel.spacing.x = unit(0.1, "lines"))
    ggsave(paste0(snakemake@params[1],"/", lin, "_cnv.png"), plot_cnv , width = 16, height = 9, dpi = 300)
    print("Plotting repeats")
    print(head(get(paste0(lin, "_repeats"))))
    plot_repeats <- ggplot(get(paste0(lin, "_repeats"))) +
        ylim(1,1)+
        geom_segment(aes(x=Start, xend=End, y = 1, yend = 1 , color = Repeat_type), linewidth = 5)+
        scale_color_manual(values = r_colors)+ 
        facet_grid(~Chromosome, scales = "free")+
        theme(panel.background = element_blank(),
            panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 0.5),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            strip.text.y = element_text(angle = 0),
            panel.spacing.x = unit(0.1, "lines"))
    print("Joining plots")
    plot <- plot_repeats / plot_cnv
    print("Saving plot")
    ggsave(paste0(snakemake@params[1],"/", lin, "_plot.png"), plot , width = 16, height = 9, dpi = 300)
}

file.create(paste0(snakemake@params[1], "/plots.done"))