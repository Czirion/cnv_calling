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
chrom_names$Lineage <- as.factor(chrom_names$Lineage)
levels(chrom_names$Lineage)
metadata <- metadata %>% 
    select(sample, Lineage = group)

print("Joining files")
cnv <- left_join(cnv, metadata, by = c("Sample" = "sample"))

cnv <- cnv %>%
    filter(Repeat_fraction == 0)

repeats <- left_join(repeats, chrom_names, by = "Accession")
print("Converting to factors")
print("Converting cnv")
cnv$Lineage <- as.factor(cnv$Lineage)
print("Converting repeats")
repeats$Lineage <- as.factor(repeats$Lineage)
repeats$Repeat_type <- ifelse(repeats$Repeat_type == "Simple_repeat", "Simple repeat", "Others")
repeats$Repeat_type <- factor(repeats$Repeat_type, levels = c("Simple repeat", "Others"))
r_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(repeats$Repeat_type))
print("Separating lineages")
lineages <- levels(cnv$Lineage)
for (lineage in lineages){
    df_cnv <- cnv %>% filter(Lineage == lineage)
    df_repeats <- repeats %>% filter(Lineage == lineage)
    cnv_name <- paste0(lineage, "_cnv")
    repeats_name <- paste0(lineage, "_repeats")
    assign(cnv_name, df_cnv)
    assign(repeats_name, df_repeats)
}
print("Done separating lineages")


set2 <- rev(brewer.pal(8, "Dark2")[1:6])
s_colors <- set2[1:length(unique(cnv$Structure))]

print("Plotting")
for (lineage in lineages){
    print(lineage)
    cnv_name <- paste0(lineage, "_cnv")
    repeats_name <- paste0(lineage, "_repeats")
    print("Plotting CNV")
    plot_cnv <- ggplot(get(cnv_name)) +
        ylim(1,1)+
        geom_segment(aes(x=Start, xend=End, y = 1, yend = 1 , color = Structure), linewidth = 5)+
        scale_color_manual(values = s_colors)+
        facet_grid(Sample~Accession, scales = "free")+
        scale_x_continuous(breaks = seq(0, max(get(cnv_name)$End, na.rm = TRUE), by = 500000), labels = scales::label_comma())+
        theme(panel.background = element_blank(),
            panel.border = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.title.x = element_blank(),
            strip.text.y = element_text(angle = 0),
            # strip.text.x = element_blank(),
            panel.spacing.x = unit(0.1, "lines"),
            panel.spacing.y = unit(0.1, "lines"))
    # print("Plotting repeats")
    # plot_repeats <- ggplot(get(repeats_name)) +
    #     ylim(1,1)+
    #     geom_segment(aes(x=Start, xend=End, y = 1, yend = 1 , color = Repeat_type), linewidth = 5)+
    #     scale_color_manual(values = r_colors)+ 
    #     facet_grid(~Chromosome, scales = "free")+
    #     theme(panel.background = element_blank(),
    #         panel.border = element_rect(colour = "lightgray", fill=NA, linewidth = 0.5),
    #         axis.title.y = element_blank(),
    #         axis.ticks.y = element_blank(),
    #         axis.text.y = element_blank(),
    #         axis.text.x = element_blank(),
    #         axis.ticks.x = element_blank(),
    #         axis.title.x = element_blank(),
    #         strip.text.y = element_text(angle = 0),
    #         panel.spacing.x = unit(0.1, "lines"))
    # print("Joining plots")
    # plot <- plot_repeats / plot_cnv
    # print("Saving plot")

    ggsave(paste0(snakemake@params[1],"/", lineage, "_cnv.png"), plot_cnv , width = 16, height = 16, dpi = 300)
    # ggsave(paste0(snakemake@params[1],"/", lineage, "_plot.png"), plot , width = 16, height = 16, dpi = 300)
}
print("Done plotting")
print("Creating plots.done file")
file.create(paste0(snakemake@params[1], "/plots.done"))