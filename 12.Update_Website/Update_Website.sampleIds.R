#### Script to update the webiste:::
####by jcb nunez
####

library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)

#### first download the current sample metadata form
samps <- fread("/Users/jcbnunez/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_25Feb2023.csv")

samps.old.dat = filter(samps, !is.na(sampleId_orig) 
                       & set != "dgn"
                       )
####

file.targ = c( "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/_includes/bed.gz.html",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/_includes/masked.sync.gz.html",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/_includes/mel.bam.html",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/_includes/pipeline_output.html",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/_includes/SNAPE.monomorphic.masked.sync.gz.html",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/bed.gz.ALL.txt ",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/masked.sync.gz.ALL.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/masked.sync.gz.GZONLY.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/masked.sync.gz.TBIONLY.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/SNAPE.monomorphic.masked.sync.gz.ALL.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/SNAPE.monomorphic.masked.sync.gz.GZONLY.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/SNAPE.monomorphic.masked.sync.gz.TBIONLY.txt",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_bed.gz.ALL.txt.csv",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_masked.sync.gz.ALL.txt.csv",
               "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_SNAPE.monomorphic.masked.sync.gz.ALL.txt.csv "
              )

###file.targ = c( 
###"/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_bed.gz.ALL.txt.csv",
###"/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_masked.sync.gz.ALL.txt.csv",
###"/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_SNAPE.monomorphic.masked.sync.gz.ALL.txt.csv "
###)

foreach(input.txt = file.targ)%do%{ #### open i
foreach(i = 1:dim(samps.old.dat)[1])%do%{ #### open i
  
old.text=samps.old.dat$sampleId_orig[i]
new.text=samps.old.dat$sampleId[i]
set=samps.old.dat$set[i]
#input.txt=file.targ

message(paste("changing",old.text, "to",new.text, "set=",  set))

if(set != "dgn"){ ### set
  
system(
 paste( "sed -i ''",
  paste("'s/", old.text, "/", new.text, "/g'", sep = "") ,
  input.txt,
  sep = " ")
)
  
#system(
#  paste( "sed -i ''",
#         paste("'s|", 
#               paste("/",old.text,"/",old.text, ".", sep = "") , 
#               "|", 
#               paste("/",new.text,"/",new.text, ".",sep = "") , 
#               "|g'", sep = ""),
#         input.txt,
#         sep = " ")
#)

} ### set
#if(set == "dgn"){ ### set
#  
#  system(
#    paste( "sed -i ''",
#           paste("'s/| ", old.text, "./| ", new.text, "./g'", sep = "") ,
#           input.txt,
#           sep = " ")
#  )
#  
#  system(
#    paste( "sed -i ''",
#           paste("'s|", 
#                 paste("DGN/",old.text,"/",old.text, ".", sep = "") , 
#                 "|", 
#                 paste("DGN/",new.text,"/",new.text, ".",sep = "") , 
#                 "|g'", sep = ""),
#           input.txt,
#           sep = " ")
#  )
#  
#} ### set



} #### close i
} #### close file targ


####

file.targ.dgn = c( 
  "/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/12.Update_Website/DGN.spec.case.html",
  "/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/12.Update_Website/DGN.download.txt",
  "/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/12.Update_Website/DGN.sync.gz.GZONLY.txt",
  "/Users/jcbnunez/Documents/GitHub/DESTv2_data_paper/12.Update_Website/DNG.sync.gz.TBIONLY.txt",
  "/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_masked.sync.gz.justDGN.csv"
)

#file.targ.dgn = c("/Users/jcbnunez/Documents/GitHub/DEST-bio.github.io/assets/util/MD5/MD5_masked.sync.gz.justDGN.csv")

samps.old.dat.dgn = filter(samps, !is.na(sampleId_orig) 
                       & set == "dgn"
)

foreach(input.txt = file.targ.dgn)%do%{ #### open i
  foreach(i = 1:dim(samps.old.dat.dgn)[1])%do%{ #### open i
    
    old.text=samps.old.dat.dgn$sampleId_orig[i]
    new.text=samps.old.dat.dgn$sampleId[i]
    set=samps.old.dat.dgn$set[i]
    #input.txt=file.targ
    
    message(paste("changing",old.text, "to",new.text, "set=",  set))
    
    #if(set != "dgn"){ ### set
    #  
    #  system(
    #    paste( "sed -i ''",
    #           paste("'s/", old.text, "/", new.text, "/g'", sep = "") ,
    #           input.txt,
    #           sep = " ")
    #  )
    #  
    #  #system(
      #  paste( "sed -i ''",
      #         paste("'s|", 
      #               paste("/",old.text,"/",old.text, ".", sep = "") , 
      #               "|", 
      #               paste("/",new.text,"/",new.text, ".",sep = "") , 
      #               "|g'", sep = ""),
      #         input.txt,
      #         sep = " ")
      #)
      
   # } ### set
    #if(set == "dgn"){ ### set
      
      system(
        paste( "sed -i ''",
               paste("'s/| ", old.text, "./| ", new.text, "./g'", sep = "") ,
               input.txt,
               sep = " ")
      )
      
      system(
        paste( "sed -i ''",
               paste("'s|", 
                     paste("DGN/",old.text,"/",old.text, ".", sep = "") , 
                     "|", 
                     paste("DGN/",new.text,"/",new.text, ".",sep = "") , 
                     "|g'", sep = ""),
               input.txt,
               sep = " ")
      )
      
    #} ### set
    
    
    
  } #### close i
} #### close file targ



