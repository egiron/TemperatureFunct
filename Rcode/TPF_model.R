#******************************************************************************
#
# Calculating grain yield using the Trapezoidal Temperature Function (TPF)
# 
# version: 1.0
# Copyright: (c) March 2022 - CIMMYT
# Authors: Azam Lashkari
#          Urs christoph schulthess
#          Ernesto Giron (e.giron.e@gmail.com)
#
#
# Last updated: December 21, 2022 - by egiron
# - Include modified Colorado NDVI and added Spain data by Azam
#
# This source is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# This code is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# A copy of the GNU General Public License is available on the World Wide Web
# at <http://www.gnu.org/copyleft/gpl.html>. You can also obtain it by writing
# to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
# MA 02111-1307, USA.
#
#******************************************************************************
rm(list=ls())
options(warn=-1) # Suppress summarise info
#options(scipen = 0, digits = 3)
# ---------------------------
# LOADING LIBRARIES
# ---------------------------
library (tidyverse, warn.conflicts = FALSE)
library(ggpubr, warn.conflicts = FALSE)
library (dplyr, warn.conflicts = FALSE)
library (zoo, warn.conflicts = FALSE)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

options(dplyr.summarise.inform = FALSE)

# ========================================================
# SETUP MODEL PARAMETERS
# ========================================================
# DATASETS
Weather_csv="Meteo_data_Final_20221219.csv" #"Meteo_STc_data_Final_20220112.csv"
NDVI_csv="NDVI_Ave_Final_20221219.csv" #"NDVI_Ave_Final_20220112.csv"
Pheno_csv="Pheno_Ave_Final_20221219.csv" #"Pheno_Ave_Final_20220112.csv"
Metrics_csv="combinations_statsxmodel_TPF.csv"
Combinations_csv="combinations_Yield_TPF.csv"

# GLOBAL VARIABLES
OUTPUT_DIR <- "results/TPF"
WD <- NULL
WeatherFile <- NULL
NDVIFile <- NULL
PhenoFile <- NULL
metricsFile <- NULL
combinationsFile <- NULL
RUE <- 3
DRYMATTER <- 0.8
FACTOR_TON_HA <- 0.01
YIELD_FACTOR <- DRYMATTER * FACTOR_TON_HA
NDVI_lowerThreshold <- 0.16
TminFactor <- 0.25
Topt <- 18
Toptmin <- 15
Toptmax <- 25
Tmin <- 9
Tmax <- 34
is_VPDStress <- TRUE

if (is_VPDStress==TRUE){
  Metrics_csv <- "combinations_statsxmodel_TPF_vpdStres.csv"
  Combinations_csv <- "combinations_Yield_TPF_vpdStress.csv"
  OUTPUT_DIR <- "results/TPF_vpdStress"
}else{
  Metrics_csv <- "combinations_statsxmodel_TPF.csv"
  Combinations_csv <- "combinations_Yield_TPF.csv"
  OUTPUT_DIR <- "results/TPF_noStress"
}

# CHARTS
EQ_FONTSIZE_CHART <- 3.5
LEGEND_FONTSIZE_CHART <- 2.6
LEGEND_POSITION_CHART_X <- 7.1
LEGEND_POSITION_CHART_Y <- 1.6
CHART_WIDTH <- 25.7
CHART_HEIGHT <- 16.1
CHART_DPI <- 320

# ==============================================================
# DO NOT MODIFY THESE FUNCTIONS IF YOU ARE NOT SURE!
# ---------------------------------
# FUNCTIONS TO SUPPORT THE MODEL
# ---------------------------------

# Function to load the datasets and create a metric files to save results.
load_datasets <-function(wd='.', Weather_csv=Weather_csv, NDVI_csv=NDVI_csv,
                         Pheno_csv=Pheno_csv, Metrics_csv=Metrics_csv,
                         Combinations_csv=Combinations_csv) {
  print(paste0("Loading datasets from folder:",wd))
  WD <<- wd
  setwd(WD)
  # Weather data
  WeatherFile <<- read.csv (paste(WD,"data/",Weather_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
  # NDVI data
  NDVIFile <<- read.csv (paste(WD,"data/",NDVI_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
  # Phenology data 
  PhenoFile <<- read.csv (paste(WD,"data/",Pheno_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
  # Create a global File to save all metrics
  if (!dir.exists("results")){
    dir.create(file.path(WD, "results"), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create directory to save all combinations in one file
  if (!dir.exists(OUTPUT_DIR)){
    dir.create(file.path(WD, OUTPUT_DIR), recursive = TRUE, showWarnings = FALSE)
  }
  # Create file to save all combinations results of crop yield (merge all of the CSVs)
  combinationsFile <<- paste(WD,OUTPUT_DIR,"/",Combinations_csv, sep="", collapse = NULL)
  
  if (file.exists(combinationsFile)) {
    unlink(combinationsFile, force=T)
  }
  # Create new one
  if (is_VPDStress==TRUE){
    df_combinationsFile <- data.frame("country"=character(), "location"=character(), "loc_code"=character(),
                                      "cycle"=numeric(), "lat"=numeric(), "lon"=numeric(), "SimYield"=numeric(), "ObsYield"=numeric(),
                                      "RUE"=numeric(), "Tmin"=numeric(), "Toptmin"=numeric(), "Toptmax"=numeric(), "Tmax"=numeric(), "TminFactor"=numeric(), "TmaxFactor"=numeric(),
                                      "Lvpd"=numeric(), "Uvpd"=numeric(), "SFvpd_Lthres"=numeric(), "SFvpd_Uthres"=numeric(),
                                      stringsAsFactors=FALSE)
  } else {
    df_combinationsFile <- data.frame("country"=character(), "location"=character(), "loc_code"=character(),
                                      "cycle"=numeric(), "lat"=numeric(), "lon"=numeric(), "SimYield"=numeric(), "ObsYield"=numeric(),
                                      "RUE"=numeric(), "Tmin"=numeric(), "Toptmin"=numeric(), "Toptmax"=numeric(), "Tmax"=numeric(), "TminFactor"=numeric(), "TmaxFactor"=numeric(),
                                      stringsAsFactors=FALSE)
  }
  write.table(df_combinationsFile, file = combinationsFile, sep = ',',
              col.names = TRUE, append = FALSE,
              row.names = FALSE, na="NA")
  
}

# Function to pre-processing the raw datasets and estimate some useful variables
preprocess_datasets <- function(NDVI_lowerThreshold, NDVIFile, PhenoFile, 
                                WeatherFile, TminFactor=0.25){
  data <- NULL
  # print("Processing NDVI")
  # First we need to interpolate NDVI and produce daily NDVI
  # remove Year, Month, Day, DOY from data frame
  NDVIFile = subset(NDVIFile, select = -c(Year, Month, Day))
  # fill missing dates, or populate missing dates
  NDVIFile <- NDVIFile %>% 
    group_by(location, loc_code, cycle) %>%
    mutate (Date = as.Date (phenotype_date)) %>%
    complete (Date = seq.Date (min (Date), max (Date), by = "day"))
  # add DOY to NDVIFile
  #DOY <- yday (NDVIFile$Date)
  DOY <- as.numeric(format(as.POSIXct(NDVIFile$Date, format = "%Y-%m-%d"), format="%j"))
  NDVIFile <- cbind(NDVIFile, DOY)
  data <- tibble(NDVIFile)
  # change column name
  data <- data %>% rename(DOY = "...10" )
  # interpolate missing NDVIs 
  inter <- data %>%
    group_by(location, loc_code, cycle) %>%
    mutate_at(vars (starts_with("NDVI")),
              funs(zoo::na.locf(zoo::na.locf(na.approx(., na.rm = FALSE, rule = 1), na.rm = FALSE),
                                fromLast = TRUE)))
  
  # write.table (inter, file = "inter.txt", sep = "\t")
  # -----------------------------------
  # print("Processing Phenology...")
  # for yield simulation we need NDVI from Heading to Maturity, 
  # so we need to get interpolated NDVI from heading to Maturity from PhenoFile
  # In PhenoFile convert sowing date to DATE format
  PhenoFile$SowingDateQC <- as.Date (PhenoFile$SowingDateQC)
  # convert days to heading and days to maturity to date (DATE), sowing date = 0
  PhenoFile = mutate (PhenoFile, 
                      Heading_date = as.Date (SowingDateQC + Days_To_Heading ),
                      Maturity_date = as.Date (SowingDateQC + Days_To_Maturity))
  # delete unwanted columns: 
  Pheno = subset(PhenoFile, select = -c( Days_To_Heading, Days_To_Maturity, ObsYield))
  # reshape from wide format to long format use gather 
  Pheno_HeadToMat <- Pheno %>%
    pivot_longer(
      cols = c(starts_with('Heading_date'), starts_with('Maturity_date')),
      names_to = "phenology", values_to = "days_date"
    )
  # fill missing dates, or populate missing dates in Pheno_HeadToMat
  Pheno_HeadToMat <- Pheno_HeadToMat %>% 
    group_by(location, loc_code, cycle) %>%
    mutate (Date = as.Date (days_date)) %>%
    complete (Date = seq.Date (min (Date), max (Date), by = "day"))
  
  # add Phenology date to NDVI data (join Pheno_HeadToMat and NDVI file) 
  pheno_ndvi <- left_join(Pheno_HeadToMat, inter, by = c("location","loc_code", "cycle", "Date"))
  # add a column to data: DAH: Days after Heading
  pheno_ndvi <- pheno_ndvi %>% group_by(location, loc_code, cycle) %>%
    mutate ( DAH = sequence(rle(cycle)$lengths) )
  
  # get last row of each group : if NDVI at maturity is NA , set the NDVI at last row NDVI = 0.16  
  pheno_ndviQC <- pheno_ndvi %>% group_by(location, loc_code, cycle) %>%
    mutate (NDVI = ifelse(row_number()==n() , NDVI_lowerThreshold, NDVI))
  
  # repeat the interpolation and fill the missing values (NAs) 
  NdviFinal <- pheno_ndviQC %>%
    group_by(location, loc_code, cycle) %>%
    mutate_at(vars (starts_with("NDVI")),
              funs(na.locf(na.locf(na.approx(., na.rm = FALSE, rule = 1), na.rm = FALSE),
                           fromLast = TRUE)))
  print("Calculating iPAR from NDVI")
  # calculate iPAR_from_ndvi
  iPAR <- mutate (NdviFinal, iPAR_fromNDVI = NDVI*1.25-0.19)
  # NB: 0.01 <= iPAR_fromNDVI <= 0.95 , so we need to replace negative values with 0.01
  iPAR["iPAR_fromNDVI"][iPAR["iPAR_fromNDVI"] < 0.01] <- 0.01
  iPAR["iPAR_fromNDVI"][iPAR["iPAR_fromNDVI"] > 0.95] <- 0.95
  
  # ------------------------------
  print("Processing Weather...")
  # add iPAR_fromNDVI to Weather data for each GID and each Environment
  # WeatherFile
  # add Date column to WeatherFile
  WeatherFile$Date <- as.Date(with(WeatherFile,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")
  # convert Date in iPAR to Date format
  iPAR$Date <- as.Date (iPAR$Date, format = "%m/%d/%Y")
  # get weather data for each location, cycle, gid (different gid has different days to heading and days to maturity)
  data <- left_join(iPAR, WeatherFile, by = c ("location", "Date"))
  # Rename the columns
  data <- rename(data,
                 SolRad = `Shortwave Radiation [MJ/m2/d]` , 
                 TMAX = `TemperatureMax [C]`,
                 TMIN = `TemperatureMin [C]`,
                 #TC = `Canopy Temperature [C]`, # The new weather file has not this variable
                 PCP = `Precipitation [mm/d]`,
                 iPAR = iPAR_fromNDVI,
                 VPDMAX = `Vapor Pressure Deficit max [kPa]`,
                 #VPDTC = `Vapor Pressure Deficit at Tc [kPa]`, # The new weather file has not this variable
                 RhumMin = `Relative Humidity min [%]`,
                 RhumMax = `Relative Humidity max [%]`,
                 WINDSPEED = `Wind Speed 2m [m/s]`)
  
  TMIN_PERC_FACTOR <- TminFactor
  TMAX_PERC_FACTOR <- round(1.0 - TminFactor, 2)
  # Estimate Tday and TC
  data <- dplyr::mutate (data, Tdaymax = TMIN_PERC_FACTOR*TMIN + TMAX_PERC_FACTOR*TMAX )
  
  return(data)
}
# ---------------------------------
# Add some stats per model
addStatstoYieldResults <-function(data){
  # number of days with TMAX > 34 C and number of days with TMIN < 9 C
  data_stats <- data %>%
    group_by(location, loc_code, cycle) %>%
    mutate (
      ndays_tmin9 = ifelse(TMIN<9, 1, 0) ,
      ndays_tmax34 = ifelse(TMAX>34, 1, 0)
    ) %>%
    summarise(ndays_tmin9 = sum(ndays_tmin9),
              ndays_tmax34 = sum(ndays_tmax34)
    )
  # Average of NDVI, iPAR, etc from heading to maturity
  data_stats2 <- data %>%
    group_by(location, loc_code, cycle) %>%
    summarise(ave_Tdaymax = mean(Tdaymax),
              ave_NDVI = mean(NDVI),
              ave_iPAR = mean(iPAR)
    )
  
  data_stats <- left_join(data_stats, data_stats2, by = c("location","loc_code", "cycle"))
  
  return(data_stats)
}

# ---------------------------------------------
# Calculate Trapezoidal Temperature Function
calculateTPF <- function(Tday, Tmin, Toptmin, Toptmax, Tmax){
  tpf <- 0
  if (Toptmin > Toptmax){
    warning("Min Optimum Temperature greater than Maximum Opt. Temperature")
    tpf <- NA
  } else if (Toptmax > Tmax){
    warning("Max Optimum Temperature greater than Maximum Temperature")
    tpf <- NA
  } else {
    if ((Tday <= Tmin) | (Tday >= Tmax)){
      tpf <- 0
    } else if ((Tday >= Tmin) & (Tday < Toptmin)){
      interceptLeft <- lsfit(x=c(Tmin,Toptmin), y=c(0,1))$coefficients[1]
      slopeLeft <- lsfit(x=c(Tmin,Toptmin), y=c(0,1))$coefficients[2]
      tpf <- Tday * slopeLeft + interceptLeft
    } else if ((Tday > Toptmax) & (Tday <= Tmax)){
      interceptRight <- lsfit(x=c(Toptmax, Tmax), y=c(1,0))$coefficients[1]
      slopeRight <- lsfit(x=c(Toptmax, Tmax), y=c(1,0))$coefficients[2]
      tpf <- Tday * slopeRight + interceptRight
    } else if ((Tday >= Toptmin) & (Tday <= Toptmax)) { 
      tpf <- 1
    }
  }
  return(tpf)
}
calculateTPF_wrap <- function( Tday, Tmin, Toptmin, Toptmax, Tmax) {
  y <- vector("numeric", length (Tday)) 
  for (i in 1:length(Tday)) {
    y[[i]] <- calculateTPF( Tday=Tday[[i]], Tmin, Toptmin, Toptmax, Tmax) 
  }
  return(y)
}

# GPP calculation
getGPP <- function(SolRad, RUE, Tm, iPAR, stressFactor=1, is_VPDStress=F){
  gpp <- vector("numeric", length (SolRad)) 
  for (i in 1:length(SolRad)) {
    if (is_VPDStress==TRUE){
      gpp[[i]] <- SolRad[[i]]*0.5*RUE*Tm[[i]]*iPAR[[i]]*stressFactor[[i]]
    } else{
      gpp[[i]] <- SolRad[[i]]*0.5*RUE*Tm[[i]]*iPAR[[i]]*stressFactor
    }
  }
  return(gpp)
}

# No Stress condition
calculateGPP <- function(df, RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax){
  # Weights for TMAX and TMIN
  TMIN_PERC_FACTOR <- TminFactor
  TMAX_PERC_FACTOR <- round(1.0 - TminFactor, 2)
  # Estimate Tday
  df <- dplyr::mutate (df, 
                       Tdaymax = TMIN_PERC_FACTOR*TMIN + TMAX_PERC_FACTOR*TMAX,
                       TPFTMAX = calculateTPF_wrap(Tdaymax, Tmin, Toptmin, Toptmax, Tmax),
                       GPPTMAX_TPF = getGPP(SolRad,RUE,TPFTMAX,iPAR)
  )
  return(df)
}

# GPP calculation VPD Stress
calculateGPP_VPDStress <- function(df, RUE, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres){
  df <- dplyr::mutate (df,
                       SFvpd = calculateSFvpd(VPDMAX, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres),
                       GPPSFvpd_TPF = getGPP(SolRad,RUE,TPFTMAX,iPAR, SFvpd, is_VPDStress=T)
  )
  
  return(df)
}

#----------------------------------------------------------
# SFvpd
calculateSFvpd <- function(VPDx, Lvpd, Uvpd, SFvpd_Lthres, SFvpd_Uthres){
  sfvpd <- vector("numeric", length (VPDx)) 
  for (i in 1:length(VPDx)) {
    if( VPDx [[i]] <= 0 ) {
      sfvpd[[i]] <- SFvpd_Lthres 
    } else if (VPDx [[i]] > 0 & VPDx [[i]] <= Lvpd) {
      sfvpd[[i]] <- SFvpd_Uthres #1   
    } else if (VPDx [[i]] > Lvpd & VPDx [[i]] < Uvpd) {
      sfvpd[[i]] <- 1 - (VPDx[[i]]-1)/(4.6-1) 
    } else {
      sfvpd[[i]] <- SFvpd_Lthres 
    }
  }
  return(sfvpd) 
}

#----------------------------------------------------------
# Yield Calculation
getYield <- function(df, phenoVars, is_VPDStress=F){
  #print("Estimating Yield...")
  df <- df %>% group_by(location, loc_code, cycle)
  if (is_VPDStress==TRUE){
    yield <- df %>% summarise( SGPPSFvpd_TPF = sum(GPPSFvpd_TPF) )
    yield <- dplyr::mutate (yield, 
                            SimYield = round(SGPPSFvpd_TPF * YIELD_FACTOR, 2)
    )
  } else {
    yield <- df %>% summarise( SGPPTMAX_TPF = sum(GPPTMAX_TPF) )
    yield <- dplyr::mutate (yield, 
                            SimYield = round(SGPPTMAX_TPF * YIELD_FACTOR, 2) 
    )
  }
  return(yield)
}

# --------------------------------------------------------
# Report & Chart
reportYield <- function(yield, RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, phenoVars, data_stats, 
                        Lvpd=NA, Uvpd=NA, SFvpd_Lthres=NA, SFvpd_Uthres=NA,
                        shrtReport=F, saveReport=F, mergeReport=T, saveFig=F, is_VPDStress=F, grp='country'){
  TFname = 'TPF'
  # Put back some variables
  yield <- left_join(yield, phenoVars, by = c("location","loc_code", "cycle"))
  yield$days_heading_to_maturity <- difftime(as.Date(yield$Days_To_Maturity), as.Date (yield$Days_To_Heading) , units = c("days"))
  yield$ObsYield <- round(yield$ObsYield, 2) 
  
  # Short Report
  if (shrtReport==TRUE){
    col_order <- c("country", "location", "loc_code", "cycle", 
                   "lat", "lon", "SimYield", "ObsYield")
    yield <- yield[, col_order]
  } else {
    yield <- left_join(yield, data_stats, by = c("location","loc_code", "cycle"))
  }
  
  # Add treatment attributes
  TxFactor = round(1.0 - TminFactor, 2)
  if (is_VPDStress==TRUE){
    yield <- mutate (yield, RUE=RUE, Tmin=Tmin, Toptmin=Toptmin, Toptmax=Toptmax, Tmax=Tmax, TminFactor=TminFactor, TmaxFactor=TxFactor, Lvpd=Lvpd, Uvpd=Uvpd, SFvpd_Lthres=SFvpd_Lthres, SFvpd_Uthres=SFvpd_Uthres)
  } else{
    yield <- mutate (yield, RUE=RUE, Tmin=Tmin, Toptmin=Toptmin, Toptmax=Toptmax, Tmax=Tmax, TminFactor=TminFactor, TmaxFactor=TxFactor)
  }
  
  if (saveReport==TRUE){
    if (is_VPDStress==TRUE){
      fname = paste(c("Yield_RUE",RUE,"_L",Lvpd,"_U",Uvpd,"_SFvpd",SFvpd_Lthres,"-",SFvpd_Uthres,"_TmnFact",TminFactor,"_TmxFact",TxFactor,"_Tm",Tmin,"_Tx",Tmax,"_Toptmin",Toptmin,"_Toptmax",Toptmax,"_",TFname,".csv"), collapse = "")
    } else {
      fname = paste(c("Yield_RUE",RUE,"_TmnFact",TminFactor,"_TmxFact",TxFactor,"_Tm",Tmin,"_Tx",Tmax,"_Toptmin",Toptmin,"_Toptmax",Toptmax,"_",TFname,".csv"), collapse = "")
    }
    
    print(paste0("Saving results to ", fname))
    write.csv (yield, file = paste(WD,"results/",fname, sep="", collapse = NULL)) 
  } 
  
  # Save charts
  if (saveFig==TRUE){
    if (is_VPDStress==TRUE){
      subtitle = paste(c("RUE: ",RUE," - VPDmx: [L:",Lvpd,", U:", Uvpd,"] - SFvpd thresholds: [",SFvpd_Lthres,", ",SFvpd_Uthres,"] - Tm: ",Tmin," - Tx:",Tmax," - Toptmin: ",Toptmin," - Toptmax: ",Toptmax," - Tmfact: ",TminFactor ), collapse = "")
      chart_fname = paste(c("Yield_RUE",RUE,"_L",Lvpd,"_U",Uvpd,"_SFvpd",SFvpd_Lthres,"-",SFvpd_Uthres,"_Tm",Tmin,"_Tx" ,Tmax,"_Toptmin",Toptmin,"_Toptmax",Toptmax,"_Tmfact",TminFactor,"_",TFname), collapse = "")
    } else {
      subtitle = paste(c("RUE: ",RUE," - Tm: ",Tmin," - Tx:" ,Tmax," - Toptmin: ",Toptmin," - Toptmax: ",Toptmax," - Tmfact: ",TminFactor), collapse = "")
      chart_fname = paste(c("Yield_RUE",RUE,"_Tm",Tmin,"_Tx" ,Tmax,"_Toptmin",Toptmin,"_Toptmax",Toptmax,"_Tmfact",TminFactor,"_",TFname), collapse = "")
    }
    chartYieldxmodel(yield, TFname, grp, subtitle, fname=chart_fname, saveFig=saveFig)
  }
  
  # Save combinations in a global file
  if (mergeReport==TRUE){
    write.table(yield, file = combinationsFile,
                sep = ',', append = T,
                col.names = !file.exists(combinationsFile),
                row.names = FALSE, na="NA")
  }
  
  return(yield)
  
}


# Calculate metrics 
calculateMetrics <- function(obs, pred){
  m_mae <- MAE(obs,pred)
  m_mse <- MSE(obs,pred)
  m_rmse <- RMSE(obs,pred)
  m_rmsre <- RMSRE(obs,pred)
  m_EF <- EF(obs,pred)
  m_MAPE <- MAPE(obs,pred)
  fit1 = lm(obs~pred) #, data = df, na.action = na.omit #Create the linear regression
  if (!is.na(summary(fit1)$r.squared) & !is.null(summary(fit1)$r.squared)){
    m_R2 <- signif(summary(fit1)$r.squared, 8) #R2(obs,pred)
  } else {
    m_R2 <- NA
  }
  if (!is.na(summary(fit1)$adj.r.squared) & !is.null(summary(fit1)$adj.r.squared)){
    m_AdjR2 <- signif(summary(fit1)$adj.r.squared, 8)
  } else {
    m_AdjR2 <- NA
  }
  if (!is.na(fit1$coef[[1]]) & !is.null(fit1$coef[[1]])){
    m_intercept <- signif(fit1$coef[[1]],5 )
  } else {
    m_intercept <- NA
  }
  if (!is.na(fit1$coef[[2]]) & !is.null(fit1$coef[[2]])){
    m_slope <- signif(fit1$coef[[2]], 5)
  } else {
    m_slope <- NA
  }
  if (!is.null(summary(fit1)$coef) && !is.na(summary(fit1)$coef)
      && length(summary(fit1)$coef) >= 8 ){
    m_pvalue <- signif(summary(fit1)$coef[2,4], 5)
  } else {
    m_pvalue <- NA
  }
  m_accuracy <- (100 - m_MAPE)
  # CCC
  cccrval <- CCC(obs,pred)
  Cb <- round(cccrval$C.b, digits = 2)
  m_CCC <- round(cccrval$rho.c[,1], digits = 2)
  # cat("MAE:", m_mae, "\n", "MSE:", m_mse, "\n",
  #     "RMSE:", m_rmse, "\n", "RMSRE:", m_rmsre, "\n",
  #     "MAPE:", m_MAPE, "\n", "CCC:", m_CCC, "\n", "Cb:", Cb, "\n",
  #     "p-value:", m_pvalue, "\n", "R-squared:", m_R2, "\n",
  #     "Adj-R-squared:", m_AdjR2, "\n", "EF:", m_EF, "\n",
  #     "Accuracy: ", m_accuracy)
  nd = 3 # Number of decimals
  return(list(round(m_mae,nd), round(m_mse,nd), round(m_rmse,nd), round(m_rmsre,nd), 
              round(m_MAPE,nd), m_pvalue, round(m_R2,nd), round(m_AdjR2,nd), round(m_EF,nd), 
              round(m_accuracy,nd), round(m_intercept,nd), round(m_slope,nd), m_CCC, Cb
  ))
}
# Create Charts per Model
# --------------------------------------
chart4Yield <- function(yield, metrics, TFname, pred_fld=pred_fld, 
                        grp='country', subtitle, fname){
  m_mae <- round(metrics[[1]],3)
  m_mse <- round(metrics[[2]],3)
  m_rmse <- round(metrics[[3]],3)
  m_rmsre <- round(metrics[[4]],3)
  m_MAPE <- round(metrics[[5]],3)
  m_pvalue <- signif(metrics[[6]], 5)
  m_R2 <- round(metrics[[7]],3)
  m_AdjR2 <- round(metrics[[8]],3)
  m_EF <- round(metrics[[9]],3)
  m_accuracy <- round(metrics[[10]],3)
  m_intercept <- round(metrics[[11]],3)
  m_slope <- round(metrics[[12]],3)
  m_CCC <- round(metrics[[13]],2)
  Cb <- round(metrics[[14]],2)
  yield <- yield[is.finite(yield[[pred_fld]]) && is.finite(yield$ObsYield), ]
  vrname <- "" #sub(paste("_",TFname, sep="", collapse = NULL),"",sub("Yield","", pred_fld))
  p <- ggplot(yield) +
    geom_smooth(data=yield,  aes(x=yield[[pred_fld]], y= ObsYield),method="lm", formula = y~x, se=F, col="red", size=0.8, linetype = "dashed", show.legend = F) +
    geom_point(aes(x=yield[[pred_fld]], y= ObsYield, color=yield[[grp]], group=yield[[grp]]), alpha=0.8,size=2) +
    labs(title = vrname,
         x="Observed Yield   ("~t ~ha^-1~")", 
         y=bquote(paste("Simulated Yield ", .(vrname), "     (", t ~ha^-1, ")"))
    ) +
    #y=paste("Simulated Yield -",vrname,u, sep="", collapse = NULL) ) +
    geom_label(aes(x = LEGEND_POSITION_CHART_X, 
                   y = LEGEND_POSITION_CHART_Y), 
               hjust = 0, 
               size=LEGEND_FONTSIZE_CHART,
               label = paste("RMSE =",m_rmse, 
                             "\nMAE =",m_mae, 
                             "\nRMSRE:", m_rmsre,
                             "\nR² = ", m_R2, 
                             "\nEF =",m_EF,
                             "\np-value =",m_pvalue,
                             "\nCCC =",m_CCC,
                             "\nAccuracy =",m_accuracy
               ))
  if (!is.na(m_intercept) && !is.na(m_slope)){
    p <- p + geom_label(
      aes(x = 0.5, y = 9.5), hjust = 0, col='black',
      label = lm_eqn3(m_intercept, m_slope), size=EQ_FONTSIZE_CHART, label.size = NA,
      parse = TRUE, show.legend=F)
  }
  p <- p + geom_abline(intercept = 0, slope = 1, colour='#d3d3d3', linetype = "solid", size = 0.5) +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
    scale_colour_discrete(name = grp) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

getMetricsAndChartYield <- function(df, obs_fld='ObsYield', TFname, pred_fld,
                                    grp, subtitle, fname){
  p_arr <- list()
  for (i in 1:length(pred_fld)) {
    metrics <- calculateMetrics(df[[obs_fld]], df[[pred_fld[i] ]])
    p <- chart4Yield(df, metrics,TFname, pred_fld=pred_fld[i], grp, subtitle, fname )
    p_arr <- c(p_arr, list(p))
  }
  return(p_arr)
}
# Create a Figure for Observed vs Simulated Yield
chartYieldxmodel <- function(df, TFname='TPF', grp='country', subtitle, fname, saveFig=F){
  pred_fld <- c("SimYield")
  p_arr <- getMetricsAndChartYield(df, obs_fld='ObsYield', TFname, pred_fld, 
                                   grp, subtitle, fname)
  p1 <- p_arr[[1]]
  #p2 <- p_arr[[2]]
  #
  p <- ggarrange(p1 + theme(plot.title = element_blank()), # + rremove( "xlab") + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank() ), 
                 #p2 + theme(plot.title = element_blank()), #theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), plot.margin = margin(l = 2, r=1) ), # + rremove("ylab") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
                 labels = NULL,
                 #labels = c("a)", "b)"),
                 ncol = 1, nrow = 1,
                 #widths = c(1,1), heights = c(1,1),  
                 align = "hv",
                 legend = "bottom", common.legend = TRUE)
  
  #plt_title <- "Observed vs Simulated Yield"
  plt_title <- paste("Observed vs Simulated Yield \n(", TFname," no stress)", sep="", collapse = NULL)
  p <- annotate_figure(p, top = text_grob(plt_title, face = "bold", size = 14),
                       bottom = text_grob(subtitle, color = "red", face = "bold", size = 10)
                       #fig.lab = "Figure 1", fig.lab.face = "bold"
  )
  # Save chart
  if (saveFig==TRUE){
    if (!dir.exists(OUTPUT_DIR)){
      dir.create(file.path(WD, OUTPUT_DIR), recursive = TRUE, showWarnings = FALSE)
    }
    ggsave(paste(WD,OUTPUT_DIR,"/",fname, ".png", sep="", collapse = NULL), 
           width = CHART_WIDTH, height = CHART_HEIGHT, units = "cm", 
           dpi=CHART_DPI, bg='white')
  }
  return(p)
}

lm_eqn3 <- function(intercept, slope){
  eq <- substitute(italic(y) == a + b %.% italic(x), 
                   list(a = intercept, b = slope ))
  as.character(as.expression(eq));
}

# ==============================================
# METRICS
# ==============================================
#
# ----------------------------------------------
# Mean absolute error (MAE)
# ----------------------------------------------
MAE <- function(obs,pred){
  mae <- mean(abs(obs - pred))
  return (mae)
}
# ----------------------------------------------
# Mean Squared Error (MSE)
# ----------------------------------------------
MSE <- function(obs,pred){
  mse <- mean((obs - pred)^2)
  return (mse)
}
# ----------------------------------------------
# Root Mean Square Error (RMSE) 
# ----------------------------------------------
RMSE <- function(obs,pred){
  rmse <- sqrt(mean((obs - pred)^2))
  return (rmse)
}
# ----------------------------------------------
# Mean Absolute Percentage Error (MAPE)
# ----------------------------------------------
MAPE = function(obs,pred){
  mape <- mean(abs((obs - pred)/obs))*100
  return(mape)
}
# ----------------------------------------------
# R-Squared
# ----------------------------------------------
R2 <- function(obs,pred){
  #r2 <- 1 - (sum((pred - obs)^2) / sum((obs - mean(obs))^2))
  r2 <- cor(obs,pred)^2
  return (r2)
}
# ----------------------------------------------
# Root Mean Square Relative Error (RMSRE)
# ----------------------------------------------
RMSRE <- function(obs,pred){
  rmsre <- 100 * sqrt(mean(((obs - pred) / obs)^2))
  return (rmsre)
}
# ----------------------------------------------
# Nash–Sutcliffe model efficiency (EF)
# ----------------------------------------------
# EF is a distance measure that compares model MSE with the MSE of using 
# the average of measured values as an estimator. Therefore, EF is useful 
# for making statements about the skill of a model relative to this simple 
# reference estimator. For a model that simulates perfectly,
# EF = 1, and for a model that has the same squared error of simulation as 
# the mean of the measurements, EF = 0. EF is positive for a model that has 
# a smaller squared error than the mean of the measurements.
EF <- function(obs,pred){
  ef <- 1 - ( sum((obs - pred)^2) / sum((obs - mean(obs))^2) )
  return(ef)
}

# ------------------------------------------------------------------
# Lin’s Concordance Correlation Coefficient (CCC)
# Computes Lin's (1989, 2000) concordance correlation coefficient for 
# agreement on a continuous measure obtained by two methods. The 
# concordance correlation coefficient combines measures of both precision 
# and accuracy to determine how far the observed data deviate from the 
# line of perfect concordance (that is, the line at 45 degrees on a square 
# scatter plot). 
# 
# Extracted from epiR: Tools for the Analysis of Epidemiological Data
# install.packages("epiR")
# https://rdrr.io/cran/epiR/
# https://search.r-project.org/CRAN/refmans/DescTools/html/CCC.html

# rho.c : the concordance correlation coefficient.
# s.shift : the scale shift.
# l.shift : the location shift.
# C.b : a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees. No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989, page 258).
# blalt : a data frame with two columns: mean the mean of each pair of measurements, delta vector y minus vector x.
# sblalt : a data frame listing the average difference between the two sets of measurements, the standard deviation of the difference between the two sets of measurements and the lower and upper confidence limits of the difference between the two sets of measurements. If rep.measure == TRUE the confidence interval of the difference is adjusted to account for repeated observations across individual subjects.
# nmissing : a count of the number of measurement pairs ignored due to missingness.

# $rho.c
# est lower upper
# 1   1     1     1
# 
# $s.shift
# [1] 1
# 
# $l.shift
# [1] 0
# 
# $C.b
# [1] 1
# 
# $blalt
# mean delta
# 1     1     0
# 2     2     0
# 3     3     0
# 4     4     0
# 5     5     0
# 6     6     0
# 7     7     0
# 8     8     0
# 9     9     0
# 10   10     0
# 
# $sblalt
# est delta.sd lower upper
# 1   0        0     0     0
# 
# $nmissing
# [1] 0

# ------------------------------------------------------------------
CCC = function(x, y, ci = "z-transform", conf.level = 0.95, rep.measure = FALSE, subjectid){
  
  N. <- 1 - ((1 - conf.level) / 2)
  zv <- qnorm(N., mean = 0, sd = 1)
  
  dat <- data.frame(x, y)
  id <- complete.cases(dat)
  nmissing <- sum(!complete.cases(dat))
  dat <- dat[id,]
  
  k <- length(dat$y)
  yb <- mean(dat$y)
  sy2 <- var(dat$y) * (k - 1) / k
  sd1 <- sd(dat$y)
  
  xb <- mean(dat$x)
  sx2 <- var(dat$x) * (k - 1) / k
  sd2 <- sd(dat$x)
  
  r <- cor(dat$x, dat$y)
  sl <- r * sd1 / sd2
  
  sxy <- r * sqrt(sx2 * sy2)
  p <- 2 * sxy / (sx2 + sy2 + (yb - xb)^2)
  
  delta <- (dat$x - dat$y)
  rmean <- apply(dat, MARGIN = 1, FUN = mean)
  blalt <- data.frame(mean = rmean, delta)
  
  # Scale shift:
  v <- sd1 / sd2
  # Location shift relative to the scale:
  u <- (yb - xb) / ((sx2 * sy2)^0.25)
  # Variable C.b is a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees (a measure of accuracy). 
  # No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989 page 258).
  # C.b <- (((v + 1) / (v + u^2)) / 2)^-1
  
  # The following taken from the Stata code for function "concord" (changed 290408):
  C.b <- p / r
  
  # Variance, test, and CI for asymptotic normal approximation (per Lin [March 2000] Biometrics 56:325-5):
  sep <- sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2)) / (r)^2 + (2 * (p)^3 * (1 - p) * (u)^2 / r) - 0.5 * (p)^4 * (u)^4 / (r)^2 ) / (k - 2))
  ll <- p - (zv * sep)
  ul <- p + (zv * sep)
  
  # Statistic, variance, test, and CI for inverse hyperbolic tangent transform to improve asymptotic normality:
  t <- log((1 + p) / (1 - p)) / 2
  set = sep / (1 - ((p)^2))
  llt = t - (zv * set)
  ult = t + (zv * set)
  llt = (exp(2 * llt) - 1) / (exp(2 * llt) + 1)
  ult = (exp(2 * ult) - 1) / (exp(2 * ult) + 1)
  
  # Calculate delta.sd if repeated measures:
  if(rep.measure == TRUE){  
    # Make sure subject is a factor:
    dat$sub <- subjectid
    if(!is.factor(dat$sub)) dat$sub <- as.factor(dat$sub)
    
    # Number of subjects:
    nsub <- length(levels(dat$sub))      
    
    # One way analysis of variance:
    model <- aov(delta ~ dat$sub)           
    
    # Degrees of freedom:
    MSB <- anova(model)[[3]][1]       
    
    # Sums of squares:
    MSW <- anova(model)[[3]][2]       
    
    # Calculate number of complete pairs for each subject:
    pairs <- NULL
    
    for(i in 1:nsub){
      pairs[i] <- sum(is.na(delta[dat$sub == levels(dat$sub)[i]]) == FALSE)
    }
    
    sig.dl <- (MSB - MSW) / ((sum(pairs)^2 - sum(pairs^2)) / ((nsub - 1) * sum(pairs)))
    delta.sd <- sqrt(sig.dl + MSW)
  }
  
  # Calculate delta.sd if no repeated measures:
  if(rep.measure == FALSE){
    delta.sd <- sqrt(var(delta, na.rm = TRUE))
  }
  
  # Upper and lower bounds for Bland Altmann plot:
  ba.p <- mean(delta)
  ba.l <- ba.p - (zv * delta.sd)
  ba.u <- ba.p + (zv * delta.sd)
  sblalt <- data.frame("est" = ba.p, "delta.sd" = delta.sd, "lower" = ba.l, "upper" = ba.u) 
  
  if(ci == "asymptotic"){
    rho.c <- data.frame(p, ll, ul)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)  
  }
  
  else if(ci == "z-transform"){
    rho.c <- data.frame(p, llt, ult)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)
  }
  
  return(rval)
}


# Get all combinations for TPF (no stress)
getAllCombinations_TPF_noStress <- function(comb_args, rtrnList=T){
  df_cmb <- do.call(expand.grid, comb_args)
  #df_cmb <-df_cmb[order(df_cmb$Tmin, df_cmb$Toptmin), ]
  # Remove combinations with Tmin > Toptmin
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Tmin < df_cmb$Toptmin)
  # Remove combinations with Toptmin > Toptmax
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Toptmin < df_cmb$Toptmax)
  # Remove combinations with Toptmax > Tmax
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Toptmax < df_cmb$Tmax)
  # Convert to List of list of inputs
  if (rtrnList==T)
    df_cmb <- split(df_cmb, seq(nrow(df_cmb)))
  return(df_cmb)
}

getAllCombinations_TPF_vpdStress <- function(comb_args, rtrnList=F){
  df_cmb <- do.call(expand.grid, comb_args)
  # Remove combinations with Lvpd >= Uvpd
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Uvpd > df_cmb$Lvpd)
  # Remove combinations with Tmin > Toptmin
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Tmin < df_cmb$Toptmin)
  # Remove combinations with Toptmin > Toptmax
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Toptmin < df_cmb$Toptmax)
  # Remove combinations with Toptmax > Tmax
  df_cmb <- df_cmb %>% dplyr::filter(df_cmb$Toptmax < df_cmb$Tmax)
  
  # df_cmb_final <-df_cmb_final[order(df_cmb_final$Tmin, df_cmb_final$Toptmin), ]
  # Convert to List of list of inputs
  if (rtrnList==T)
    df_cmb <- split(df_cmb, seq(nrow(df_cmb)))
  return(df_cmb)
} 


# Create Dataframe with statistics per model
createDF <- function(yield, RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax, pred_fld, is_VPDStress=F,
                     Lvpd=NA, Uvpd=NA, SFvpd_Lthres=NA, SFvpd_Uthres=NA
){
  TFname = "TPF"
  metrics <- calculateMetrics(yield$ObsYield, yield[[pred_fld]])
  TmaxFactor <- round(1.0 - TminFactor, 2)
  nd <- 3
  if (is_VPDStress==TRUE){
    df_stats <- data.frame(
      RUE = c(RUE), 
      TminFactor = c(TminFactor),
      TmaxFactor = c(TmaxFactor),
      Tmin = c(Tmin), 
      Toptmin = c(Toptmin),
      Toptmax = c(Toptmax),
      Tmax = c(Tmax), 
      Lvpd = c(Lvpd), 
      Uvpd = c(Uvpd), 
      SFvpd_L = c(SFvpd_Lthres), 
      SFvpd_U = c(SFvpd_Uthres), 
      TempFunct = c(TFname),
      MAE = c(round(metrics[[1]],nd)),
      MSE = c(round(metrics[[2]],nd)),
      RMSE = c(round(metrics[[3]],nd)),
      RMSRE = c(round(metrics[[4]],nd)),
      MAPE = c(round(metrics[[5]],nd)),
      pvalue = c(signif(metrics[[6]], 5)),
      R2 = c(round(metrics[[7]],2)),
      AdjR2 = c(round(metrics[[8]],2)),
      EF = c(round(metrics[[9]],nd)),
      intercept = c(round(metrics[[11]],nd)),
      slope = c(round(metrics[[12]],nd)),
      Cb = c(round(metrics[[14]],2)),
      CCC = c(round(metrics[[13]],2)),
      Accuracy = c(round(metrics[[10]],2)),
      stringsAsFactors=FALSE)
  } else {
    df_stats <- data.frame(
      RUE = c(RUE), 
      TminFactor = c(TminFactor),
      TmaxFactor = c(TmaxFactor),
      Tmin = c(Tmin), 
      Toptmin = c(Toptmin),
      Toptmax = c(Toptmax),
      Tmax = c(Tmax), 
      TempFunct = c(TFname),
      MAE = c(round(metrics[[1]],nd)),
      MSE = c(round(metrics[[2]],nd)),
      RMSE = c(round(metrics[[3]],nd)),
      RMSRE = c(round(metrics[[4]],nd)),
      MAPE = c(round(metrics[[5]],nd)),
      pvalue = c(signif(metrics[[6]], 5)),
      R2 = c(round(metrics[[7]],2)),
      AdjR2 = c(round(metrics[[8]],2)),
      EF = c(round(metrics[[9]],nd)),
      intercept = c(round(metrics[[11]],nd)),
      slope = c(round(metrics[[12]],nd)),
      Cb = c(round(metrics[[14]],2)),
      CCC = c(round(metrics[[13]],2)),
      Accuracy = c(round(metrics[[10]],2)),
      stringsAsFactors=FALSE)
  }
  
  return(df_stats)
}

saveMetricsResults <- function(df, fname){
  if (!dir.exists(OUTPUT_DIR)){
    dir.create(file.path(WD, OUTPUT_DIR), recursive = TRUE, showWarnings = FALSE)
  }
  m = paste(WD, OUTPUT_DIR, "/",fname, sep="", collapse = NULL)
  if (file.exists(m)) {
    unlink(m, force=T)
  }
  write.table(df, file = m, sep = ';', append = T,
              col.names = !file.exists(m), 
              row.names = FALSE, na="NA")
}

# Function to run simulation and get metrics in parallel
metrics_combinations_TPF_noStress <- function(inputs, data, phenoVars, data_stats, saveReport=F, mergeReport=T, saveFig=F){
  RUE = inputs[['RUE']]
  Tmin = inputs[['Tmin']]
  Toptmin = inputs[['Toptmin']]
  Toptmax = inputs[['Toptmax']]
  Tmax = inputs[['Tmax']]
  TminFactor = inputs[['TminFactor']]
  TFname <- "TPF"
  pred_fld <- c("SimYield")
  dataTPF <- calculateGPP(data, RUE=RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax)
  yieldTPF <- getYield(dataTPF, phenoVars)
  r_yieldTPF <- reportYield(yieldTPF, RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, phenoVars, data_stats, 
                             Lvpd=NA, Uvpd=NA, SFvpd_Lthres=NA, SFvpd_Uthres=NA, shrtReport=T, 
                             saveReport=saveReport, mergeReport=mergeReport, saveFig=saveFig)
  df_stats <- createDF(r_yieldTPF, RUE=RUE, TminFactor=TminFactor, Tmin=Tmin, Toptmin=Toptmin, 
                       Toptmax=Toptmax, Tmax=Tmax, pred_fld=pred_fld)
  # Save metrics in a global file
  # write.table(df_stats, file = metricsFile, 
  #             sep = ';', append = T,
  #             col.names = !file.exists(metricsFile), 
  #             row.names = FALSE, na="NA")
  
  return(df_stats)
}

metrics_combinations_TPF_vpdStress <- function(inputs, data, phenoVars, data_stats, saveReport=F, mergeReport=T, saveFig=F){
  RUE = inputs[['RUE']]
  Tmin = inputs[['Tmin']]
  Toptmin = inputs[['Toptmin']]
  Toptmax = inputs[['Toptmax']]
  Tmax = inputs[['Tmax']]
  TminFactor = inputs[['TminFactor']]
  Lvpd = inputs[['Lvpd']]
  Uvpd = inputs[['Uvpd']]
  SFvpd_Lthres = inputs[['SFvpd_Lthres']]
  SFvpd_Uthres = inputs[['SFvpd_Uthres']]
  TFname <- "TPF"
  pred_fld <- c("SimYield")
  dataTPF <- calculateGPP(data, RUE=RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax)
  dataTPF <- calculateGPP_VPDStress(dataTPF, RUE=RUE, Lvpd=Lvpd, Uvpd=Uvpd, 
                                     SFvpd_Lthres=SFvpd_Lthres, SFvpd_Uthres=SFvpd_Uthres)
  yieldTPF <- getYield(dataTPF, phenoVars, is_VPDStress=T)
  
  r_yieldTPF <- reportYield(yieldTPF, RUE=RUE, Tmin=Tmin, Toptmin=Toptmin, Toptmax=Toptmax, Tmax=Tmax, 
                            TminFactor=TminFactor, phenoVars, data_stats, 
                            Lvpd=Lvpd, Uvpd=Uvpd, SFvpd_Lthres=SFvpd_Lthres, SFvpd_Uthres=SFvpd_Uthres, shrtReport=T, 
                            saveReport=saveReport, mergeReport=mergeReport, saveFig=saveFig, is_VPDStress=T)
  
  df_stats <- createDF(r_yieldTPF, RUE=RUE, TminFactor=TminFactor, Tmin=Tmin, Toptmin=Toptmin, 
                       Toptmax=Toptmax, Tmax=Tmax, pred_fld=pred_fld, is_VPDStress=T, 
                       Lvpd=Lvpd, Uvpd=Uvpd, SFvpd_Lthres=SFvpd_Lthres, SFvpd_Uthres=SFvpd_Uthres)
  
  # Save metrics in a global file
  # write.table(df_stats, file = metricsFile, 
  #             sep = ';', append = T,
  #             col.names = !file.exists(metricsFile), 
  #             row.names = FALSE, na="NA")
  
  return(df_stats)
}

# ====================== END FUNCTIONS =======================


# =============================================================
# =============================================================
# RUN MODEL 
# =============================================================

# ---------------------------
# Setup working directory
# ---------------------------
#wd <- getwd() 
wd <- "~/Desktop/CIMMYT/AzamLashkari_2021/VDP and temperature functions/models/"
setwd(wd)
# ---------------------------
# LOADING DATASETS
# ---------------------------
load_datasets(wd, Weather_csv, NDVI_csv, Pheno_csv, Metrics_csv, Combinations_csv)
phenoVars <- subset(PhenoFile, select = c( location, loc_code, cycle, country, Days_To_Heading, Days_To_Maturity, ObsYield, lat, lon))
# ---------------------------
# PROCESSING DATA
# ---------------------------
data <- preprocess_datasets(NDVI_lowerThreshold, NDVIFile, PhenoFile, WeatherFile, TminFactor)
data_stats <- addStatstoYieldResults(data)

# --------------------------------
# TESTING NONE STRESS CONDITIONS
# --------------------------------
# Estimate GPP
# dataTPF <- calculateGPP(data, RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax)
# # Estimate Yield
# yieldTPF <- getYield(dataTPF, phenoVars)
# # Report & Chart
# r_yieldTPF <- reportYield(yieldTPF, RUE, Tmin, Toptmin, Toptmax, Tmax, TminFactor, phenoVars, data_stats, shrtReport=T, saveReport=T, mergeReport=T, saveFig=T)
# # Display Figure
# chartYieldxmodel(r_yieldTPF, TFname="TPF", grp='loc_code', subtitle="TPF", fname="TPF", saveFig=F)
# ------------- END TESTING -------------------

# ==============================================================
# Run combinations in parallel
# ==============================================================
# Combinations for RUE, TminFactor, Tmin, Toptmin, Toptmax, Tmax in TPF (no stress)
combination_args_TPF_noStress <- list(
  RUE = c(3.0), #c(seq(2.8, 3.2, 0.1)),
  Tmin = c(0:10),
  Toptmin = c(10:25),
  Toptmax = c(25:30),
  Tmax = c(30:45),
  TminFactor = c(0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5)
)

# Test
# combination_args_TPF_noStress <- list(
#   RUE = c(3.0),
#   Tmin = c(9.0),
#   Toptmin = c(15.0),
#   Toptmax = c(25.0),
#   Tmax = c(34.0),
#   TminFactor = c(0.25)
# )

# Get all combinations for TPF
cmb_list_TPF_noStress <- getAllCombinations_TPF_noStress(combination_args_TPF_noStress)

start <- proc.time()
cl <- parallel::makeCluster( parallel::detectCores() - 1,  type = "FORK" ) #"PSOCK" #FORK
doParallel::registerDoParallel(cl = cl)
df_result_combinations_TPF_noStress <- foreach( i = cmb_list_TPF_noStress, .packages = c('dplyr'), 
                                                .combine = 'rbind' ) %dopar% {
  tryCatch({
    metrics_combinations_TPF_noStress(i, data, phenoVars, data_stats, saveReport=F, mergeReport=T)
  }, error = function(e) return(paste0("Error occurred for '", i, "'", " The error is '", e, "'")))
}
saveMetricsResults(df_result_combinations_TPF_noStress, "combinations_TPF_noStress.csv")
parallel::stopCluster(cl = cl)
print(proc.time() - start)
# user    system   elapsed 
# 19918.102  6271.394  179032.481   -> ~50 hours -> ~2 days
# ========================== END RUN IN PARALLEL ====================================
# combinations_TPF_noStress <- read.csv (paste(WD,OUTPUT_DIR,"/",Combinations_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
# View(combinations_TPF_noStress)



# 2) Combinations for Tmin, Toptmin, Toptmax, Tmax in TPF (VPD stress)
combination_args_TPF_vpdStress <- list(
  RUE = c(3.0), #c(seq(2.8, 3.2, 0.1)),
  Tmin = c(0:10),
  Toptmin = c(10:25),
  Toptmax = c(15:30),
  Tmax = c(30:45),
  TminFactor = c(0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5),
  Lvpd = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5), 
  Uvpd = c(1, 1.5, 2, 2.5, 3, 3.5, 4),
  SFvpd_Lthres = c(0.2, 0.4, 0.6, 0.8), 
  SFvpd_Uthres = c(1)
)

# Test
combination_args_TPF_vpdStress <- list(
  RUE = c(3.0),
  Tmin = c(9.0),
  Toptmin = c(15.0),
  Toptmax = c(25.0),
  Tmax = c(34.0),
  TminFactor = c(0.25),
  Lvpd = c(1.0), 
  Uvpd = c(4.0),
  SFvpd_Lthres = c(0.2), 
  SFvpd_Uthres = c(1)
)

cmb_list_TPF_vpdStress <- getAllCombinations_TPF_vpdStress(combination_args_TPF_vpdStress, rtrnList=F) # Get all combinations

start <- proc.time()
cl <- parallel::makeCluster( parallel::detectCores() - 1, type = "FORK" ) #"PSOCK" #FORK -> avoid duplicated values or send data several time to GPU 
doParallel::registerDoParallel(cl = cl)
# Using the dataframe
df_result_combinations_TPF_vpdStress <- foreach( i = 1:nrow(cmb_list_TPF_vpdStress), .packages = c('dplyr'), .combine = 'rbind' ) %dopar% {
  tryCatch({
    r <- cmb_list_TPF_vpdStress[i, ]
    metrics_combinations_TPF_vpdStress(r, data, phenoVars, data_stats, saveReport=F, mergeReport=T)
  }, error = function(e) return(paste0("Error occurred for '", i, "'",
                                       " The error is '", e, "'")))
}
# df_result_combinations_TPF_vpdStress
saveMetricsResults(df_result_combinations_TPF_vpdStress, "combinations_TPF_vpdStress.csv")
parallel::stopCluster(cl = cl)
print(proc.time() - start)
# user  system elapsed 
# 
# For 7620480 obs of 8 variables
# (((617.445 * 7620480) / 258720) / 60) / 24 -> ~13 days in PC 4 cores

# combinations_TPF_vpdStress <- read.csv (paste(WD,OUTPUT_DIR,"/",Combinations_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
# View(combinations_TPF_vpdStress)



# -------------------------
# READ MERGED RESULTS
# -------------------------
# This is not possible in some systems, due to RAM memory. 
# It will produce: Error: vector memory exhausted (limit reached?)
# combinations_Yield_TPF <- read.csv (paste(WD,OUTPUT_DIR,"/",Combinations_csv, sep="", collapse = NULL), check.names = FALSE, stringsAsFactors = FALSE)
# View(combinations_Yield_TPF)




