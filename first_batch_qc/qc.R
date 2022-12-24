### libraries
  library(data.table)
  library(patchwork)
  library(ggplot2)
  library(googlesheets4)

### mapping data
  setwd("/Users/alanbergland/Documents/GitHub/DESTv2_data_paper")
  dat <- fread("DESTv2_data_paper/first_batch_qc/Mapping.stats.final.txt")
  dat[,c("Sample", "Group", "coverage_across_reference_Dmel"), with=F][order(Group)]
  setnames(dat, "Sample", "sampleId")

### QC data
  gs <- "https://docs.google.com/spreadsheets/d/1bf6x_ow7KqUSNQAWOLkzJhwI86Spw65eXwl5dY-sWag/edit#gid=152737368"
  qual <- as.data.table(read_sheet(ss=gs))
  #qual <- qual[2:253,]
  setnames(qual, names(qual), c("Library_Name", "conc", "addvol", "finalvol", "tot", "DIN",
        "Result_Original", "Status", "libtype", "conc2", "concNm", "size", "result", "comment", "totalreads", "q30"))
  setnames(qual, "Library_Name", "SequencingId")

  qual.dt <- data.table(SequencingId=unlist(qual$SequencingId[2:254]),
                        DIN=unlist(qual$DIN[2:254]),
                        Result_Original=unlist(qual$Result_Original[2:254]),
                        Status=unlist(qual$Status[2:254]),
                        result=unlist(qual$result[2:254]),
                        totalreads=unlist(qual$totalreads[2:254]))

### metadata
  samps <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2_samples.csv")

### merges
  qual.dt <- merge(qual.dt, samps, by="SequencingId", all=T)[set=="DrosEU_3"]
  qual.dt <- merge(qual.dt, dat, by="sampleId", all.x=T)
  str(qual.dt)

  qual.dt[,totalreads:=as.numeric(as.character(totalreads))]

  qual.dt[,exp_coverage:=totalreads*(1-PERCENT_DUPLICATION) * 200/(180000000)]

### save
  save(qual.dt, file="/Users/alanbergland/Documents/GitHub/DESTv2_data_paper/first_batch_qc/qc.Rdata")

### basic pllot
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads, color=Result_Original)) + geom_point()
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads, color=Status)) + geom_point()
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads, color=result)) + geom_point() + facet_grid(~Result_Original)

  ggplot(data=qual.dt, aes(y=coverage_across_reference_Holo, x=totalreads, color=result)) + geom_point() + facet_grid(~Result_Original)

  ggplot(data=qual.dt, aes(x=coverage_across_reference_Holo, y=coverage_across_reference_Dmel, color=result)) + geom_point() + facet_grid(~Result_Original)

  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads*(1-PERCENT_DUPLICATION), color=result)) +
  geom_point()

### accross mel
  tot_read_plot_nolabel <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads, color=result)) +
  geom_point()


  exp_read_plot_nolabel <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads*(1-PERCENT_DUPLICATION), color=result)) +
  geom_point()

  tot_read_plot <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads, color=result)) +
  geom_point() +
  geom_text_repel(
    data=qual.dt[coverage_across_reference_Dmel<20],
    aes(label = sampleId),
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "grey50"
  )


  exp_read_plot <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Dmel, x=totalreads*(1-PERCENT_DUPLICATION), color=result)) +
  geom_point() +
  geom_text_repel(
    data=qual.dt[coverage_across_reference_Dmel<20],
    aes(label = sampleId),
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .15,
    nudge_y = .5,
    color = "grey50"
  )


### across the hologenome
  tot_read_plot_holo <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Holo, x=totalreads, color=result)) +
  geom_point() +
  geom_text_repel(
    data=qual.dt[coverage_across_reference_Dmel<20],
    aes(label = sampleId),
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .015,
    nudge_y = .05,
    color = "grey50"
  )


  exp_read_plot_holo <-
  ggplot(data=qual.dt, aes(y=coverage_across_reference_Holo, x=totalreads*(1-PERCENT_DUPLICATION), color=result)) +
  geom_point() +
  geom_text_repel(
    data=qual.dt[coverage_across_reference_Dmel<20],
    aes(label = sampleId),
    size = 3,
    min.segment.length = 0,
    seed = 42,
    box.padding = 0.5,
    max.overlaps = Inf,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .015,
    nudge_y = .05,
    color = "grey50"
  )


  layout <- "
  AB
  CD"

  mega <-
  tot_read_plot_nolabel + exp_read_plot_nolabel +
  tot_read_plot + exp_read_plot +
  plot_annotation(tag_level="A") +
  plot_layout(guides = 'collect', design=layout)

  ggsave(mega, file="~/mega_destv2.pdf", h=10, w=10)
