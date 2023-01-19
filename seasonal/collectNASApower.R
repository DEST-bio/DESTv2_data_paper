# ijob -A berglandlab_standard -c20 -p standard --mem=40G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R


### libraries
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)

## load
  fl <- list.files("~/DESTv2_data_paper/seasonal/nasapower_sampleId/", full.names=T)

  power.dt <- foreach(fl.i=fl)%dopar%{
    # fl.i <- fl[1]
    message(fl.i)
    load(fl.i)
    return(power.dt)
  }
  power.dt <- rbindlist(power.dt)

### save
  save(power.dt, file="/scratch/aob2x/nasaPower.allpops.Rdata")
