# ijob -A berglandlab_standard -c5 -p standard --mem=20G
# module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
	library(data.table)
	library(foreach)
	library(rnoaa)
	library(sp)
  library(lubridate)
  library(doMC)
	library(climate)
  registerDoMC(1)
  library(readxl)
  library(nasapower)

### cli
	args <- commandArgs()
	message(pasge("args: ", args, sep=""))
	jobID <- as.numeric(args[1])
	# jobID <- 24
	message(jobID)

### load in samps
  setwd("~")
  #setwd("/Users/alanbergland/Documents/GitHub")
	samps <- fread("DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")
	samps <- samps[!is.na(jday)]
	samps <- samps[jobID]

	samps

### power function
	getPower <- function(i, yearOffset=0) {
			daily_single_ag <- get_power(
			  community = "ag",
			  lonlat = c(samps[i]$long, samps[i]$lat),
			  pars = c("RH2M", "T2M", "PRECTOTCORR"),
			  dates = c(paste(samps[i]$year+yearOffset, "-01-01", sep=""), paste(samps[i]$year+yearOffset, "-12-31", sep="")),
			  temporal_api = "hourly",
				time_standard="UTC"
			)
			daily_single_ag <- as.data.table(daily_single_ag)
			daily_single_ag
	}

### use NASA power
	power.dt <- list()
	counter <- 1
	for(i in c(1:dim(samps)[1])) {
    # i <- 1
		for(j in c(-1, 0, 1)) {
			message(paste(i, j, dim(samps)[1]))

			power.dt[[counter]] <- getPower(i, yearOffset=j)
			Sys.sleep(5)
			power.dt[[counter]][,sampleId:=samps[i]$sampleId]
			power.dt[[counter]][,date:=ymd_hms(paste(paste(YEAR, MO, DY, sep="-"), paste(HR, ":00:00", sep=""), sep=" "))]
			power.dt[[counter]][,queryDate:=Sys.time()]
			power.dt[[counter]][,yearOffset:=j]
			counter <- counter + 1
		}
	}

	power.dt <- rbindlist(power.dt)


	save(power.dt, file=paste("~/DESTv2_data_paper/seasonal/nasapower_sampleId/", samps$sampleId, ".Rdata", sep=""))
