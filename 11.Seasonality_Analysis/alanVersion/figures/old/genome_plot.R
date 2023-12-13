### scp aob2x@rivanna.hpc.virginia.edu:~/dest2_glm_baypass_annotation.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)

### load
  load("~/dest2_glm_baypass_annotation.Rdata")

###
  setkey(m, col)
  fisher.test(table(m[J(c("missense_variant", "synonymous_variant"))]$col,
                    m[J(c("missense_variant", "synonymous_variant"))]$C2_std_median>5.5))

  m[J("stop_gained")][C2_std_median>5.5]
  prop.table(table(m[C2_std_median>5.5]$col))

  ggplot(data=m, aes(x=af, y=C2_std_median)) + geom_hex()


### Adh
  c2 <- ggplot(data=m[chr=="2L"][pos>14e6 & pos<15e6], aes(x=pos, y=C2_std_median)) +
        geom_vline(xintercept=14617051) + geom_hline(yintercept=5.5)+ geom_line(color="grey", alpha=.5) +
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr)
  glm <- ggplot(data=m[chr=="2L"][pos>14e6 & pos<15e6], aes(x=pos, y=-log10(p_lrt))) +
        geom_vline(xintercept=14617051) + geom_line(color="grey", alpha=.5)+
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr)
  xtx <- ggplot(data=m[chr=="2L"][pos>14e6 & pos<15e6], aes(x=pos, y=XtXst_median)) +
        geom_vline(xintercept=14617051) + geom_line(color="grey", alpha=.5)+
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr)

  layout <- "
  A
  B
  C"

  c2 + glm + xtx + plot_layout(design=layout)

### genome
  c2 <- ggplot(data=m[C2_std_median>5.5], aes(x=pos, y=C2_std_median)) +
        geom_hline(yintercept=5.5)+ geom_line(color="grey", alpha=.5) +
        geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr) +
        geom_point(data=m[C2_std_median>5.5][rnp<.05][XtXst_median>200])

  glm <- ggplot(data=m[rnp<.05], aes(x=pos, y=-log10(p_lrt))) +
        geom_line(color="grey", alpha=.5)+
        geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr) +
        geom_point(data=m[C2_std_median>5.5][rnp<.05][XtXst_median>200])

  xtx <- ggplot(data=m[XtXst_median>200], aes(x=pos, y=XtXst_median)) +
        geom_line(color="grey", alpha=.5)+
        geom_vline(data=candidates, aes(xintercept=pos_mean), color="green") +
        geom_point(size=.75, alpha=.5, color="red") + facet_grid(~chr) +
        geom_point(data=m[C2_std_median>5.5][rnp<.05][XtXst_median>200])


  layout <- "
  A
  B
  C
  D"

  c2 + glm + xtx + c2.wza.plot + plot_layout(design=layout)

table(m[C2_std_median>5.5][rnp<.05][XtXst_median>200]$col)
  m[C2_std_median>5.5][rnp<.05][XtXst_median>200][col=="missense_variant"]


m[XtXst_median>400][chr=="3L"][pos>20e6]

  m[XtXst_median>250][C2_std_median>5.5][pos>14e6 & pos<15e6]
  mean(m[chr=="2L"][pos==14617051]$XtXst_median> m$XtXst_median)
  m[gene=="Adh"][which.max(XtXst_median)]
  mean(m[gene=="Adh"][which.max(XtXst_median)]$XtXst_median> m$XtXst_median)

  prod <- ggplot(data=m, aes(x=pos, y=-log10(C2.rnp*rnp))) + geom_point(size=.75, alpha=.5) + facet_grid(~chr) + geom_line()


  #glm <- ggplot(data=m[-log10(p_lrt)>3.5 & C2_std_median>5], aes(x=pos, y=-log10(p_lrt)), size=.75, alpha=.5) + geom_point() + facet_grid(~chr)

  c2 / glm
  m[chr=="2L"][which.max(C2_std_median)]
m[-log10(p_lrt)>3.5 & C2_std_median>5][which.max(C2_std_median)]
m[-log10(p_lrt)>3.5 & C2_std_median>5][which.min(p_lrt)]
m[which.min(C2.rnp*rnp)]
candidates <- data.table(genes=unique(m[C2_std_median>10]$gene))
universe <- data.table(genes=unique(m$gene))

write.table(candidates, file="~/candidates.csv", quote=F, row.names=F, col.names=F, sep=",")
write.table(universe, file="~/universe.csv", quote=F, row.names=F, col.names=F, sep=",")
