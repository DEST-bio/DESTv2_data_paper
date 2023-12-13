system("scp aob2x@rivanna.hpc.virginia.edu:~/XtX_C2_glm.windows.Rdata ~/.")

### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)

  setwd("/Users/alanbergland/Documents/GitHub")

### load window data
  load("DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures/XtX_C2_glm.windows.Rdata")

### candidate genes
  candidates <- data.table(chr=c("3R", "3L", "3L", "2L", "2L", "2R", "2L"),

                          pos_mean=c(13241247, 3317423, 20221145, 14617051, 7558557, 14875642, 18260220),
                          label=c("Ace", "Drs", "obst-F", "Adh", "Cyp4d21", "Cyp6a8", "?"))

### sliding window plot
nSNP_thr <- 150
  xtx.median.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") + geom_text(data=candidates, aes(y=220, x=pos_mean/1e6, label=label), size=2) +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=xtx_median), color="red") +
    facet_grid(~chr, scales="free_x") + ylab("median XtXst") + xlab("Pos (Mb)") +
     theme(axis.text.x = element_text(color = "grey20", size = 6))

  xtx.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-xtx.wZa.p), color="red") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(xtx wZa p)") + xlab("Pos (Mb)") +
     theme(axis.text.x = element_text(color = "grey20", size = 6))

  c2.median.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=(c2_median)), color="red") +
    facet_grid(~chr, scales="free_x") + ylab("median C2") + xlab("Pos (Mb)") +
     theme(axis.text.x = element_text(color = "grey20", size = 6))

  c2.wza.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
    geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") +
    geom_line(data=win.out,aes(x=pos_mean/1e6, y=-(C2.wZa.p)), color="red") +
    facet_grid(~chr, scales="free_x") + ylab("-log10(C2 wZa p)") + xlab("Pos (Mb)") +
     theme(axis.text.x = element_text(color = "grey20", size = 6))

    glm.rnp.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
      geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") +
      geom_line(data=win.out,aes(x=pos_mean/1e6, y=-log10(glm.rnp.binom.p)), color="red") +
      facet_grid(~chr, scales="free_x") + ylab("GLM RNP 5%") + xlab("Pos (Mb)") +
       theme(axis.text.x = element_text(color = "grey20", size = 6))

    glm.wZa.plot <- ggplot(data=win.out[nSNPs>nSNP_thr]) +
      geom_vline(data=candidates, aes(xintercept=pos_mean/1e6), color="green") +
      geom_line(data=win.out,aes(x=pos_mean/1e6, y=-(wZa.p)), color="red") +
      facet_grid(~chr, scales="free_x") + ylab("-log10(GLM wZa p)") + xlab("Pos (Mb)") +
       theme(axis.text.x = element_text(color = "grey20", size = 6))



  layout <- "
  ACE
  BDF"

  layout <- "
  AB
  CD
  EF"

  mega_window <- xtx.median.plot + xtx.wza.plot + c2.median.plot + c2.wza.plot + glm.rnp.plot  + glm.wZa.plot +
   plot_layout(design=layout) + plot_annotation(tag_level="a")

  ggsave(mega_window, file="DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures/glm_baypass.png", h=10, w=7)

  xtx.wza.plot+ geom_text(data=candidates, aes(y=3500, x=pos_mean, label=label)) + c2.wza.plot + glm.wZa.plot + plot_layout(ncol=1)



### candidate windows
  load("~/dest2_glm_baypass_annotation.Rdata")
  candidate <- win.out[chr=="2L"][which.min(C2.wZa.p)]
  candidate <- win.out[which.max(xtx_median)]
  candidate<- win.out[chr=="3L"][pos_mean<=5e6][which.min(C2.wZa.p)]

  fbgn <- as.data.table(fread("http://ftp.flybase.org/releases/FB2023_05/precomputed_files/genes/gene_map_table_fb_2023_05.tsv.gz"))
  fbgn[,chr:=tstrsplit(sequence_loc, ":")[[1]]]
  fbgn[,start:=as.numeric(tstrsplit(tstrsplit(sequence_loc, ":")[[2]], "\\.")[[1]])]
  fbgn[,end:=as.numeric(tstrsplit(tstrsplit(tstrsplit(sequence_loc, ":")[[2]], "\\.")[[3]], "\\(")[[1]])]
  setnames(fbgn, "##organism_abbreviation", "sp")
  drslike_pos<- fbgn[grepl("Drs", current_symbol)][sp=="Dmel"]


  cSNPs[,annot:=col]
  cSNPs[col%in%c("downstream_gene_variant", "upstream_gene_variant"), annot:="intergenic"]
  cSNPs[col%in%c("intron_variant", "splice_region_variant&intron_variant", "splice_region_variant", "splice_region_variant&synonymous_variant"), annot:="intronic"]
  table(cSNPs$annot)
  setkey(m, chr, pos)
  cSNPs <- m[J(data.table(chr=candidate$chr, pos=candidate$pos_min:candidate$pos_max, key="chr,pos")), nomatch=0]

  

  C2_SNPs <- ggplot(data=cSNPs, aes(x=pos, y=C2_neglogp_mean)) + geom_line() + geom_hline(yintercept=quantile(m$C2_neglogp_mean, .995)) +
          geom_point(data=cSNPs[!annot%in%c("intergenic", "intronic")], aes(color=annot))
  GLM_SNPs <- ggplot(data=cSNPs, aes(x=pos, y=-log10(p_lrt))) +  geom_line() + geom_hline(yintercept=quantile(-log10(m$p_lrt), .995)) +
        geom_point(data=cSNPs[!annot%in%c("intergenic", "intronic")], aes(color=annot))
  XTX_SNPs <- ggplot(data=cSNPs, aes(x=pos, y=XtXst_median)) +   geom_line() + geom_hline(yintercept=quantile(m$XtXst_median, .995)) +
      geom_point(data=cSNPs[!annot%in%c("intergenic", "intronic")], aes(color=annot)) +
      geom_segment(data=drslike_pos, aes(y=75, yend=75, x=start, xend=end))



  layout <- "
  A
  B
  C"

  drs_mega <- C2_SNPs + GLM_SNPs + XTX_SNPs + plot_layout(design=layout) + plot_annotation(title="Drs region")
  ggsave(drs_mega, file="DESTv2_data_paper/11.Seasonality_Analysis/alanVersion/figures//Drs_combined.pdf")
