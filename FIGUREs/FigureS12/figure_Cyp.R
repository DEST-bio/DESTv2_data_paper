library(data.table)
library(gggenes)
library(ggthemes)
library(cowplot)
library(ggplot2)
library(dplyr)

fst <- fread("RESULTS/hierFST_k8_10.csv")

fgt_all_NA_SA<-ggplot(data = fst[fst$chrom=="2R" & fst$continents=="4-North_America:4-South_America", ], aes(x = chromStart, y = FGT)) +
  geom_rect(xmin = 12162496, xmax = 12259496, ymin = -Inf, ymax = Inf,
            fill = "lightblue", alpha = 0.2) +
  geom_line() + scale_x_continuous(expand = expansion(mult = c(0, .01)), labels = scales::unit_format(unit = "M", scale = 1e-6)) + theme_base() + labs(x="", title="North America (4) vs. South America (4)", y = expression(F[GT])) +
  theme(plot.background = element_rect(
    color = "white"
  ))
fgt_all_NA_SA
fgt_all_SA_E<-ggplot(data = fst[fst$chrom=="2R" & fst$continents=="4-South_America:8-Europe", ], aes(x = chromStart, y = FGT)) +
  geom_rect(xmin = 12162496, xmax = 12259496, ymin = -Inf, ymax = Inf,
            fill = "lightblue", alpha = 0.2) +
  geom_line() + scale_x_continuous(expand = expansion(mult = c(0, .01)), labels = scales::unit_format(unit = "M", scale = 1e-6)) + theme_base() + labs(x="", title="South America (4) vs. Europe (8)", y = expression(F[GT])) +
  theme(plot.background = element_rect(
    color = "white"
  ))
fgt_all_SA_E
fgt_all_NA_EU<-ggplot(data = fst[fst$chrom=="2R" & fst$continents=="4-North_America:8-Europe", ], aes(x = chromStart, y = FGT)) +
  geom_rect(xmin = 12162496, xmax = 12259496, ymin = -Inf, ymax = Inf,
            fill = "lightblue", alpha = 0.2) +
  geom_line() + scale_x_continuous(expand = expansion(mult = c(0, .01)), labels = scales::unit_format(unit = "M", scale = 1e-6)) + theme_base() + labs(x="", title="North America (4) vs. Europe (8)", y = expression(F[GT])) +
  theme(plot.background = element_rect(
    color = "white"
  ))
fgt_all_NA_EU
fgt_all<-fgt_all_NA_EU

fst_plot<-ggplot(data = fst[fst$chrom=="2R" & fst$continents=="4-North_America:8-Europe", ], aes(x = chromStart, y = FGT)) + geom_line() + scale_x_continuous(expand = expansion(mult = c(0, .01)), labels = scales::unit_format(unit = "M", scale = 1e-6)) + theme_base() + labs(x="",  y = expression(F[GT])) + xlim(c(12162496,12259496)) + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x  = element_blank()) +
  theme(plot.background = element_rect(
    color = "white"
  ))
fst_plot

genes <- fread("genes.txt")
genes_filtered <- genes[grepl("^Cyp", genes$gene), ]
genes$kk<-genes$gene
genes$color<-"gray"
genes[!grepl("^Cyp", genes$gene), ]$kk <- ""
genes[grepl("^Cyp", genes$gene), ]$color <- "lightblue"

gene_plot<-ggplot(genes, aes(xmin = start, xmax = end, y = molecule, group = gene, fill=color, forward = orientation)) +
  geom_gene_arrow(position = position_dodge(width = 2), arrow_body_height = grid::unit(5, "mm"), arrowhead_height = grid::unit(5, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_identity() +
  geom_text(data = genes %>% mutate(start = (start + end)/2), aes(x = start, y = molecule, label = kk), position = position_dodge(width = 2), size = 3,  fontface = "italic",vjust = -2) + theme_genes() + guides(fill = F) + labs(y="", x="Chromosome 2R")+ theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(), axis.text.y = element_blank()) + theme(axis.text.x = element_text(size = 14, color = "black")) + theme(axis.title.x = element_text(face = "bold", size = 16)) + scale_x_continuous(expand = expansion(mult = c(0, .01)), labels = scales::unit_format(unit = "M", scale = 1e-6), lim = c(12162496,12259496)) +
  theme(plot.background = element_rect(
    color = "white"
  )) 
gene_plot

panelA <- plot_grid(fgt_all, fst_plot,gene_plot, align = "hv", nrow = 3, rel_heights = c(0.5,1,0.3))
ggsave(filename="FGT_EU_NA_Cyp6g1.pdf", plot_grid(fgt_all, fst_plot,gene_plot, align = "hv", nrow = 3, rel_heights = c(0.5,1,0.2)), width = 23, height = 13)

panelB<-plot_grid(fgt_all_NA_EU,fgt_all_NA_SA,fgt_all_SA_E, align = "hv", nrow = 3)
ggsave(filename="FGT_all.pdf", plot_grid(fgt_all_NA_EU,fgt_all_NA_SA,fgt_all_SA_E+labs(x="Chromosome 2R"), align = "hv", nrow = 3), width = 15, height = 6)

plot_grid(panelB,panelA,ncol=1,labels="AUTO",rel_heights = c(1,0.7))

ggsave(filename="FigureS12.pdf", width = 15, height = 22)
