
samp2 <- vroom("/Users/jcbnunez/Documents/GitHub/DESTv2/populationInfo/dest_v2_samples.csv")
qcs <- vroom("/Users/jcbnunez/Documents/GitHub/0.NoGIT_Download_DrosEUdata/5.Verify.bams/multi.qc.batch1.txt")
names(qcs)[1] = "sampleId"

bads <- c("DrosEu-240", "DrosEu-18", "DrosEu-103", "DrosEu-118", "DrosEu-211", "DrosEu-231", "DrosEu-90", "DrosEu-41", "DrosEu-88", "DrosEu-194", "DrosEu-130", "DrosEu-117", "DrosEu-129", "DrosEu-145", "DrosEu-202", "DrosEu-227", "DrosEu-199", "DrosEu-223", "DrosEu-140", "DrosEu-143", "DrosEu-138", "DrosEu-99", "DrosEu-142", "DrosEu-251", "DrosEu-248")

samp2 %>%
  mutate(BAD_BATCH = case_when(SequencingId %in% bads ~ "BAD",
                               TRUE ~ "REG")) ->
  samp2.bad

left_join(qcs, samp2.bad[c("sampleId","BAD_BATCH")]) -> qcs.batch

write.table(qcs.batch, file = "/Users/jcbnunez/Documents/GitHub/0.NoGIT_Download_DrosEUdata/5.Verify.bams/qcs.batch.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")