# ijob -A berglandlab_standard -c1 -p standard --mem=20G
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

### load in samps
  setwd("~")
	samps <- fread("DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")

### only use samps with exact collection date
  samps <- samps[exactDate==T]
  table(samps$set)

### only use localities with more than one sample per year
  samps.ag <- samps[,list(.N, daysDelta=max(jday) - min(jday)), list(locality, year)][N>1][daysDelta>60]

  setkey(samps, locality, year)
  setkey(samps.ag, locality, year)
  samps <- merge(samps, samps.ag)

  table(samps$set)

### which samples were in the core20
  drosrtec <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/elife-67577-supp1-v2.xlsx"))

  samps <- merge(samps, drosrtec[,c("Sample", "Core20")], by.x="sampleId_orig", by.y="Sample", all.x=T)
  samps[Core20=="yes", Core20:=T]
  samps[is.na(Core20) | Core20=="no", Core20:=F]

  unique(samps$Core20)

  table(samps$Core20, samps$set)

### test two weather stats
	getPower <- function(i) {
			daily_single_ag <- get_power(
			  community = "ag",
			  lonlat = c(samps[i]$long, samps[i]$lat),
			  pars = c("RH2M", "T2M", "PRECTOTCORR"),
			  dates = c(paste(samps[i]$year, "-01-01", sep=""), paste(samps[i]$year, "-12-31", sep="")),
			  temporal_api = "hourly",
				time_standard="UTC"
			)
			daily_single_ag <- as.data.table(daily_single_ag)
			daily_single_ag
	}

	getLCD <- function(i) {
		tmp <- nearest_stations_nooa(country = (samps[i]$COUNTRY),
		  date = Sys.Date(),
		  add_map = F,
		  point = c(samps[i]$long, samps[i]$lat),
		  no_of_stations = 5
		)
		tmp <- as.data.table(tmp)
		setnames(tmp, "distance [km]", "dist")
		tmp[,begin:=ymd(BEGIN)]
		tmp[,end:=ymd(END)]
		tmp <- tmp[begin<=samps[i]$Date][end >=samps[i]$Date][order(dist)]


		station <- paste(tmp$USAF, tmp$WBAN, sep="")
		message(paste(samps[i]$sampleId, samps[i]$year, i, dim(samps)[1], sep=" / "))
		lcd.dt <- foreach(j=1:length(station), .errorhandling="remove")%do%{
			#lcd.tmp <- as.data.table(lcd(station=station[j], year=samps[j]$year))

			path <- rnoaa:::lcd_get(station=station[j], year=samps[i]$year)
			lcd.tmp <- fread(path)
			setnames(lcd.tmp, names(lcd.tmp), tolower(names(lcd.tmp)))
			lcd.dt <- lcd.tmp[,c("date", "hourlydrybulbtemperature", "hourlyrelativehumidity", "hourlyprecipitation", "hourlywindspeed"), with=F]
			lcd.dt[,station:=station[j]]
			lcd.dt
		}
		lcd.dt <- rbindlist(lcd.dt)
		lcd.dt[,date:=ymd_hms(gsub("T", " ", date))]
		lcd.dt[,temp:=as.numeric(as.character(hourlydrybulbtemperature))]
		lcd.dt[,humidity:=as.numeric(as.character(hourlyrelativehumidity))]
		lcd.dt[,precip:=as.numeric(as.character(hourlyprecipitation))]
		lcd.dt[,wind:=as.numeric(as.character(hourlywindspeed))]
		lcd.dt[,sampleId:=samps[i]$sampleId]
		lcd.dt

	}

### use NASA power
	power.dt <- list()

	for(i in c(1:dim(samps)[1])) {
    # i <- 1
		message(paste(i, dim(samps)[1]))
		power.dt[[i]] <- getPower(i)
		Sys.sleep(5)
		power.dt[[i]][,sampleId:=samps[i]$sampleId]
		power.dt[[i]][,date:=ymd_hms(paste(paste(YEAR, MO, DY, sep="-"), paste(HR, ":00:00", sep=""), sep=" "))]
		power.dt[[i]][,queryDate:=Sys.time()]
	}

	power.dt <- rbindlist(power.dt)

	save(power.dt, file="~/DESTv2_data_paper/dest_v2.nasa_power_weather.raw.Rdata")

### summarize in two weeks prior to sampling
  load(file="~/DESTv2_data_paper/dest_v2.nasa_power_weather.raw.Rdata")

  power.dt.ag <- power.dt[,list(max_temp=max(T2M), min_temp=min(T2M), ave_temp=mean(T2M)), list(year, sampleId, yday=yday(ymd(paste(YEAR, MO, DY, sep="-"))))]

  m <- merge(power.dt.ag, samps, by="sampleId")
  m[,daysPrior:=jday - yday]
  m <- m[daysPrior> -14 & daysPrior<=0]
  m.ag <- m[,list(prMax=mean(max_temp>30), prMin=mean(min_temp<5), ave=mean(ave_temp)), list(sampleId)]

  samps <- merge(samps, m.ag, by="sampleId")

### save
  save(samps, file="~/DESTv2_data_paper/samps_weather_summary.Rdata")


### library
  library(data.table)
  library(ggplot2)

  load("DESTv2_data_paper/samps_weather_summary.Rdata")
  samps[,sr_season:=factor(sr_season, levels=c("spring", "fall", "winter", "frost"))]


  ggplot(data=samps[sr_season%in%c("spring", "fall")], aes(x=sr_season, y=jday, group=interaction(year, locality))) +
  geom_point() +
  geom_line() +
  facet_grid(~Core20)


  ggplot(data=samps) +
  geom_line(aes(x=jday, y=locality, group=locality)) +
  geom_point(aes(x=jday, y=locality)) +
  facet_grid(Core20~year)
