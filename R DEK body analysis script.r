# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
# Christopher Vogel                                                        Konstanz, September 2018
# 
# This script was written for the analysis of CellProfiler data from a high-throughput screen
# and is separated into four parts: 
# 1. Import and concatenation of single csv tables
# 2. Creation of summaries two reduce complexity
# 3. Quality control on position summaries of all replicates
# 4. Calculation of Z-scores to score hits
# -------------------------------------------------------------------------------------------------------------------------

# 1. Import and concatenation of single csv tables
# loading of custom packages
library(plyr)           # for manipulation of data frames
library(tidyverse)   # for manipulation of data frames
library(Hmisc)  
library(broom)       # for tidy and glance functions
library(data.table)  # for manipulation of data tables
library(mmand)     # for closing and opening morphological operations
library(dtplyr)        # for plyr-like manipulation of data tables

# Here, the working directory is set to W:/
setwd("W:/")

# in case the script should be run from windows command line:
# C:\R\R-3.4.3\bin\x64\Rscript W:\script.R

# a custom function to fill up gaps in DEK body sequences and extract the longest DEK 
# body sequence within a nucleus track
find_seq <- function(x, y) {
  closed <- closing(x, c(1, 1, 1))
  opened <- opening(closed, c(1, 1, 0))
  opened[opened > 0] <- 1
  closed <- opened * closed
  closed[closed < 1] <- NA
  if (all(is.na(closed))) {
    return(0)
  } else
    closed <- na.contiguous(closed) %>% as.vector()
    closed[closed < y] <- NA
  if (all(is.na(closed))) {
    return(0)
  } else
    return(length(na.contiguous(closed)))
}

# a custom function to fill up gaps in DEK body sequences and calculate the mean number of 
# DEK bodies of the longest DEK body-containing sequence within a nucleus track
meandek <- function(x) {
  closed <- closing(x,c(1,1,1))
  opened <- opening(closed,c(1,1,0))
  opened[opened>0] <- 1
  closed <- opened * closed
  closed[closed<1] <- NA
  if (all(is.na(closed))) {
    return(0) 
  } else  
    return(mean(na.contiguous(closed)))
}

# creates variables with metadata missing from the CellProfiler csv tables
name <- "Plate-01"
rep <- 2
plate <- paste0("May-Vogel-DEK-Bodies-", name, "-batch1-0", rep)
plate_name <- paste0(gsub("-0", "", name), "-Replicate", rep)
# creates variable with regularly used Metadata
meta <- c("Metadata_Plate", "replicate", "Metadata_BasePathChris", "Metadata_SubPathChris",
  "Metadata_siRNA", "Metadata_treatment", "Metadata_WellNr", "Metadata_Pos")

# for each image a separate filtered_nuclei.csv was created by CellProfiler 
# concatenates filtered_nuclei.csv from the same plate (must be in same directory) 
# into a single data frame
files <- list.files(path = "analysis/", pattern = "_nuclei.csv", full.names = TRUE)
temp <- lapply(files, fread, sep = ",")
fnuclei <- rbindlist(temp, fill = TRUE)
fnuclei$PlateName <- plate_name
fnuclei$replicate <- rep

# adds a new position/track identifier, that is unique for each nucleus track 
# ("PlatesiRNAWellPosTrack")
fnuclei$PlatesiRNAWellPosTrack <-
  with(fnuclei, paste0(Metadata_Plate, "--", Metadata_siRNA, "--W",  as.character(Metadata_WellNr),
      "--P", as.character(Metadata_Pos), "--TRACK", as.character(TrackObjects_Label_15)))
fnuclei %<>% select(PlatesiRNAWellPosTrack, everything()) 
fnuclei %<>% setorder(PlatesiRNAWellPosTrack, ImageNumber)

# image-based measurements are not in the filtered_nuclei.csv 
# concatenates all image.csv for image quality measurement 
# "ImageQuality_PowerLogLogSlope" for later quality control
files <- list.files(path = "analysis/", pattern = "Image.csv", full.names = TRUE)
temp <- lapply(files, fread, sep = ",")
image <- rbindlist(temp, fill = TRUE)
image$PlateName <- plate_name
image$PlatesiRNAWellPos <- with(image, paste0(Metadata_Plate, "--", Metadata_siRNA, "--W",
      as.character(Metadata_WellNr), "--P", as.character(Metadata_Pos)))
image %<>% select(PlatesiRNAWellPos, everything())
image$PlatesiRNAWellPos <- as.character(image$PlatesiRNAWellPos)

# ------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
# 2. Creation of summaries two reduce complexity
# summarizes by group siRNA/Well/Position over all time points 
# -> for each position/time lapse only one row remains
# reduces also track dimension!
image_sum <- image %>% group_by(PlatesiRNAWellPos) 
image_sum %<>% summarise(
    FileName_Tracks = unique(FileName_Tracks),
    FileName_Outlines = unique(FileName_Outlines),
    Metadata_siRNA = unique(Metadata_siRNA),
    Metadata_treatment = unique(Metadata_treatment),
    Metadata_WellNr = unique(Metadata_WellNr),
    Metadata_Pos = unique(Metadata_Pos),
    Metadata_Plate = unique(Metadata_Plate),
    PlateName = unique(PlateName),
    Metadata_BasePathChris = unique(Metadata_BasePathChris),
    Metadata_SubPathChris = unique(Metadata_SubPathChris),
    "95quant_PowerLogLog" = quantile(ImageQuality_PowerLogLogSlope_gfpdekRAW, probs = 0.95),
    upQuart_PowerLogLog = quantile(ImageQuality_PowerLogLogSlope_gfpdekRAW, probs = 0.75),
    mean_PowerLogLog = mean(ImageQuality_PowerLogLogSlope_gfpdekRAW),
    max_PowerLogLog = max(ImageQuality_PowerLogLogSlope_gfpdekRAW),
    median_PowerLogLog = median(ImageQuality_PowerLogLogSlope_gfpdekRAW),
    sum_allnuclei = sum(Count_NucleiForBG),
    sum_filtnuclei = sum(Count_FilteredNuclei)
  )

# summarise fnuclei by group siRNA/Well/Position over all time points
# -> for each position/time lapse only one row remains, reduces also track dimension!
# important measurements: 
# number of tracks for each position (num_tracks),
# sum of nuclei (sum_nuclei), mean length of tracks (mean_tracklength)
nucleisum <- fnuclei %>% group_by_at(vars = meta)
nucleisum %<>% summarise(num_tracks = max(TrackObjects_Label_15), 
  sum_nuclei = sum(ImageNumber >= 0), num_timepoints = max(Metadata_Time))
nucleisum$mean_tracklength <- nucleisum$sum_nuclei / nucleisum$num_tracks

# summarises fnuclei by position/track identifier PlatesiRNAWellPosTrack 
# over all time points/nuclei 
# -> for each track only one row remains
# important measurements:
# sum_dek: number of DEK bodies detected in each track
# mean_dekbody_count_seq: mean number of DEK bodies per DEK body sequence
# mean_mean/integ_Int_nuclei: average nuclear fluorescence intensity per track
# mean_nuclei_area: average area of detected nuclei
# tracklifetime: length of each track
# mean_intensity_dek: average fluorescence intensity of detected DEK bodies
# ratio_trackone/two/etcdek: fraction of frames within a DEK body-positive 
# track that contain at least one/two/etc DEK bodies
tracksummaries <- fnuclei %>% group_by(PlatesiRNAWellPosTrack) 
tracksummaries %<>% summarise(
    Metadata_siRNA = unique(Metadata_siRNA),
    Metadata_treatment = unique(Metadata_treatment),
    Metadata_WellNr = unique(Metadata_WellNr),
    Metadata_Pos = unique(Metadata_Pos),
    Metadata_Plate = unique(Metadata_Plate),
    trackID = unique(TrackObjects_Label_15),
    Metadata_BasePathChris = unique(Metadata_BasePathChris),
    Metadata_SubPathChris = unique(Metadata_SubPathChris),
    replicate = unique(replicate),
    tracklength_trackswdek = find_seq(Children_dekbodies_Count, 1),
    tracklength_morethanonedek = find_seq(Children_dekbodies_Count, 2),
    tracklength_morethantwodek = find_seq(Children_dekbodies_Count, 3),
    tracklength_morethanthreedek = find_seq(Children_dekbodies_Count, 4),
    tracklength_morethanfourdek = find_seq(Children_dekbodies_Count, 5),
    has_trackswdek = any(find_seq(Children_dekbodies_Count, 1)),
    has_trackswmorethanonedek = any(find_seq(Children_dekbodies_Count, 2)),
    has_trackswmorethantwodek = any(find_seq(Children_dekbodies_Count, 3)),
    has_trackswmorethanthreedek = any(find_seq(Children_dekbodies_Count, 4)),
    has_trackswmorethanfourdek = any(find_seq(Children_dekbodies_Count, 5)),
    has_dek = any(Children_dekbodies_Count),
    mean_meanInt_nuclei = mean(Intensity_MeanIntensity_gfpdek),
    mean_integInt_nuclei = mean(Intensity_IntegratedIntensity_gfpdek),
    mean_nuclei_area = mean(AreaShape_Area),
    mean_dekbody_count_seq = meandek(Children_dekbodies_Count),
    tracklifetime = n(),
    mean_intensity_dek = mean(Mean_dekbodies_Intensity_MeanIntensity_gfpdek, na.rm = TRUE)
  )
tracksummaries$ratio_trackonedek <- tracksummaries$tracklength_morethanonedek / 
  tracksummaries$tracklength_trackswdek
tracksummaries$ratio_tracktwodek <- tracksummaries$tracklength_morethantwodek /
  tracksummaries$tracklength_trackswdek
tracksummaries$ratio_trackthreedek <- tracksummaries$tracklength_morethanthreedek /
  tracksummaries$tracklength_trackswdek
tracksummaries$ratio_trackfourdek <- tracksummaries$tracklength_morethanfourdek / 
  tracksummaries$tracklength_trackswdek

# summarises track summaries by position over all tracks
# -> for each position only one row remains
# important measurements:
# sum_has_dek: number of frames that contain 
# at least one DEK body detected for each position
# sum_hasno_dek: number of frames that contain no DEK body at all
# sum_tracks: number of nuclei tracks
# fraction_of_tracks_with_dek: fraction of DEK-body positive tracks
possummaries <- tracksummaries %>% group_by_at(vars = meta)
possummaries %<>% summarise(
    mean_tracklifetime = mean(tracklifetime + 1),
    median_tracklifetime = median(tracklifetime + 1),
    mean_mean_nuclei_area = mean(mean_nuclei_area),
    mean_mean_meanInt_nuclei = mean(mean_meanInt_nuclei),
    mean_mean_integInt_nuclei = mean(mean_integInt_nuclei),
    sum_has_trackswdek = sum(has_trackswdek),
    sum_hasno_trackswdek = sum(!has_trackswdek),
    sum_tracklength_trackswdek = sum(tracklength_trackswdek),
    sum_has_dek = sum(has_dek),
    sum_hasno_dek = sum(!has_dek)
  )
possummaries$sum_tracks <- possummaries$sum_has_dek + possummaries$sum_hasno_dek
possummaries$fraction_of_tracks_with_dek <- possummaries$sum_has_trackswdek /
  possummaries$sum_tracks
possummaries$fraction_of_tracks_without_dek <- possummaries$sum_hasno_trackswdek / 
  possummaries$sum_tracks
possummaries$mean_bodyduration_tracks <- possummaries$sum_tracklength_trackswdek /
  possummaries$sum_has_trackswdek

# like possumaries, but only with tracks containing at least two consecutive 
# DEK body positive frames
deksummaries <- tracksummaries %>% filter(tracklength_trackswdek>0)
deksummaries %<>% group_by_at(vars = meta)
deksummaries %<>% summarise(
    median_tracklength_trackswdek = median(tracklength_trackswdek),
    mean_tracklength_trackswdek = mean(tracklength_trackswdek),
    mean_tracklength_morethanonedek = mean(tracklength_morethanonedek),
    mean_tracklength_morethantwodek = mean(tracklength_morethantwodek),
    mean_tracklength_morethanthreedek = mean(tracklength_morethanthreedek),
    mean_tracklength_morethanfourdek = mean(tracklength_morethanfourdek),
    mean_ratio_trackonedek = mean(ratio_trackonedek),
    mean_ratio_tracktwodek = mean(ratio_tracktwodek),
    mean_ratio_trackthreedek = mean(ratio_trackthreedek),
    mean_ratio_trackfourdek = mean(ratio_trackfourdek),
    sum_has_trackswdek = sum(has_trackswdek),
    mean_mean_dekbody_count_seq = mean(mean_dekbody_count_seq, na.rm = TRUE),
    sum_trackswmorethanonedek = sum(has_trackswmorethanonedek),
    sum_trackswmorethantwodek = sum(has_trackswmorethantwodek),
    sum_trackswmorethanthreedek = sum(has_trackswmorethanthreedek),
    sum_trackswmorethanfourdek = sum(has_trackswmorethanfourdek),
    mean_mean_nuclei_area_cellswdek = mean(mean_nuclei_area),
    mean_mean_meanInt_nuclei_cellswdek = mean(mean_meanInt_nuclei),
    mean_mean_integInt_nuclei_cellswdek = mean(mean_integInt_nuclei)
  )
deksummaries$fraction_of_trackswmorethanonedek <-
  deksummaries$sum_trackswmorethanonedek / deksummaries$sum_has_trackswdek
deksummaries$fraction_of_trackswmorethantwodek <-
  deksummaries$sum_trackswmorethantwodek / deksummaries$sum_has_trackswdek
deksummaries$fraction_of_trackswmorethanthreedek <-
  deksummaries$sum_trackswmorethanthreedek / deksummaries$sum_has_trackswdek
deksummaries$fraction_of_trackswmorethanfourdek <-
  deksummaries$sum_trackswmorethanfourdek / deksummaries$sum_has_trackswdek

# create summary statistics for each time point, keep position information
# summarise fnuclei by group siRNA/Well/Position/Time 
# -> for each time point only one row remains, reduces also track dimension!
# important measurements:
# sum_has_dek: number of frames that contain 
# at least one DEK body detected for each position
# sum_hasno_dek: number of frames that contain no DEK body at all
# sum_tracks: number of nuclei tracks
# fraction_of_tracks_with_dek: fraction of DEK-body positive tracks
nucleisummaries <- fnuclei %>% group_by(meta, Metadata_Time)
nucleisummaries %<>% summarise(
    sum_nuclei = sum(ImageNumber >= 0),
    sum_dek = sum(Children_dekbodies_Count),
    has_dek = any(Children_dekbodies_Count),
    dek_pos_nuclei = sum(Children_dekbodies_Count > 0),
    dek_neg_nuclei = sum(Children_dekbodies_Count == 0),
    nucleiwmorethanonedek = sum(Children_dekbodies_Count > 1),
    nucleiwmorethantwodek = sum(Children_dekbodies_Count > 2),
    nucleiwmorethanthreedek = sum(Children_dekbodies_Count > 3),
    nucleiwmorethanfourdek = sum(Children_dekbodies_Count > 4),
    intensity_dek = mean(Mean_dekbodies_Intensity_MeanIntensity_gfpdek, na.rm = TRUE),
    nuclei_area = mean(AreaShape_Area)
  )
nucleisummaries$fraction_of_frames_with_dek <-
  nucleisummaries$dek_pos_nuclei / nucleisummaries$sum_nuclei
nucleisummaries$mean_dekbody_count <-
  nucleisummaries$sum_dek / nucleisummaries$dek_pos_nuclei
nucleisummaries$fraction_of_nucleiwmorethanonedek <-
  nucleisummaries$nucleiwmorethanonedek / nucleisummaries$dek_pos_nuclei
nucleisummaries$fraction_of_nucleiwmorethantwodek <-
  nucleisummaries$nucleiwmorethantwodek / nucleisummaries$dek_pos_nuclei
nucleisummaries$fraction_of_nucleiwmorethanthreedek <-
  nucleisummaries$nucleiwmorethanthreedek / nucleisummaries$dek_pos_nuclei
nucleisummaries$fraction_of_nucleiwmorethanfourdek <-
  nucleisummaries$nucleiwmorethanfourdek / nucleisummaries$dek_pos_nuclei
  
nucleisummaries$PlatesiRNAWellPos <- with(nucleisummaries, paste0(Metadata_Plate, "--", 
  Metadata_siRNA, "--W", as.character(Metadata_WellNr), "--P", as.character(Metadata_Pos)))

# calculates linear regressions of cell number vs time over all positions of one well, 
# reorganizes into a data table
fitted_models <- nucleisummaries %>% group_by(PlatesiRNAWellPos)
fitted_models %<>% do(model = lm(sum_nuclei ~ Metadata_Time, data = .))
linreg <- rowwise(fitted_models) %>% tidy(model)
R2 <- rowwise(fitted_models) %>% glance(model)
linreg_wide <- dcast(setDT(linreg), PlatesiRNAWellPos ~ term, 
  value.var = c("estimate", "std.error", "statistic", "p.value"))
linreg_wide <- cbind(linreg_wide, R2["r.squared"])

# merges all above created summaries to a single summary
deksummaries_merged <- merge(deksummaries, possummaries, by = meta, all = TRUE)
deksummaries_merged <- merge(deksummaries_merged, nucleisum, by = meta, all = TRUE)
deksummaries_merged <- merge(deksummaries_merged, image_sum, by = meta, all = TRUE)
deksummaries_merged <- merge(deksummaries_merged, linreg_wide, by = "PlatesiRNAWellPos")

# adds new column RowCol with row and column number of each position for easier finding 
# on the plate layout
deksummaries_merged$Metadata_WellNr <- as.integer(deksummaries_merged$Metadata_WellNr)
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(1:12)] <-
  as.character("A")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(13:24)] <-
  as.character("B")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(25:36)] <-
  as.character("C")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(37:48)] <-
  as.character("D")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(49:60)] <-
  as.character("E")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(61:72)] <-
  as.character("F")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(73:84)] <-
  as.character("G")
deksummaries_merged$Row[deksummaries_merged$Metadata_WellNr %in% c(85:96)] <-
  as.character("H")
deksummaries_merged$Col <- (deksummaries_merged$Metadata_WellNr - 1) %% 12 + 1
deksummaries_merged$RowCol <- with(deksummaries_merged, paste0(Row, "-", Col))

# saves summarised tables, they are much easier to handle because of their reduced size
write.csv(nucleisummaries, paste0("/tables/", plate_name, "_nucleisummaries.csv"))
write.csv(deksummaries_merged, paste0("/tables/", plate_name, "_deksummaries_merged.csv"))
write.csv(tracksummaries, paste0("/tables/", plate_name, "_tracksummaries.csv"))
write.csv(image_sum, paste0("/tables/", plate_name, "_iq.csv"))
write.csv(possummaries, paste0("/tables/", plate_name, "_possummaries.csv"))

# ------------------------------------------------------------------------------------------------------------------------
# concatenates deksummaries_merged.csv from different plates (must be in same directory)
# into a single data frame
all_plates = c("Plate1-Replicate1", "Plate1-Replicate2", "Plate2-Replicate1", "Plate2-Replicate2", 
  "Plate3-Replicate1", "Plate3-Replicate4", "Plate4-Replicate1", "Plate4-Replicate2",
  "Plate5-Replicate1", "Plate5-Replicate2", "Plate5-Replicate3", "Plate6-Replicate1",
  "Plate6-Replicate2", "Plate7-Replicate2", "Plate7-Replicate3", "Plate8-Replicate1",
  "Plate8-Replicate2")
files <- list.files(path = "deksummaries/", pattern = "merged.csv", full.names = TRUE)
temp <- lapply(files, fread, sep = ",")
deksummaries_merged_QC <- rbindlist(temp, fill = TRUE)

# replaces "unknown" and wrongly annotated gene names with correct HGNC names 
# (as of October 2018)
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s227305" | Metadata_siRNA == "s56826", 
  "ATR", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s59082" | Metadata_siRNA == "s57368", 
  "ATRX", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s225845" | Metadata_siRNA=="s47345", 
  "CENPX", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s3642" | Metadata_siRNA=="s3644", 
  "CSNK2B", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA =="s37000" | Metadata_siRNA=="s37001", 
  "CTC1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s32631" | Metadata_siRNA=="s32632", 
  "DDX24", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s51989" | Metadata_siRNA=="s51990", 
  "KMT5A", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s27462" | Metadata_siRNA=="s27463", 
  "KMT5B", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s39429" | Metadata_siRNA=="s39430", 
  "KMT5C", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s35540" | Metadata_siRNA=="s35542", 
  "LBHD1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s8960" | Metadata_siRNA=="s8961", 
  "MRE11", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s30362" | Metadata_siRNA=="s30363", 
  "MRM3", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s30287" | Metadata_siRNA=="s30288", 
  "MTPAP", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s223607" | Metadata_siRNA=="s23982", 
  "MTREX", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s9892" | Metadata_siRNA=="s9893", 
  "ORC1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s24167" | Metadata_siRNA == "s24168", 
  "ORC3", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s9898" | Metadata_siRNA == "s9900", 
  "ORC4", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s9901" | Metadata_siRNA == "s9903", 
  "ORC5", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s24164" | Metadata_siRNA == "s24165", 
  "ORC6", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s18861" | Metadata_siRNA == "s18863", 
  "PCLAF", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s38873" | Metadata_siRNA == "s38874" | Metadata_siRNA == "s38875", 
  "PYM1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s20341" | Metadata_siRNA == "s20342", 
  "RACK1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s55074" | Metadata_siRNA == "s55075", 
  "RAD21L1", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s39523" | Metadata_siRNA == "s39524", 
  "RIOX2", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s9520" | Metadata_siRNA == "s9521", 
  "TONSL", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(
  Metadata_siRNA == "s37740" | Metadata_siRNA == "s37741", 
  "TRMT1L", Metadata_treatment))

deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s227428", "SRP9", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s56478", "LLPH", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s57221", "ATM", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s61393", "PRKDC", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s9896", "ORC2", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s52787", "PABPC4", Metadata_treatment))
deksummaries_merged_QC %<>% mutate(Metadata_treatment = ifelse(Metadata_siRNA ==
                                                                 "s228087", "TOP3B", Metadata_treatment))

# re-orders columns
deksummaries_merged_QC$PlateName <-
  factor(deksummaries_merged_QC$PlateName, levels = all_plates)
# deletes columns that are not interesting
deksummaries_merged_QC <- select(deksummaries_merged_QC, -c(
  Metadata_BasePathChris, Metadata_SubPathChris, sum_trackswmorethanonedek,
  sum_trackswmorethantwodek, sum_trackswmorethanthreedek,
  sum_trackswmorethanfourdek, sum_trackswmorethanfivedek,
  num_tracks, FileName_Tracks, FileName_Outlines, upQuart_PowerLogLog,
  median_PowerLogLog, Row, Col, PathName_Outlines, PathName_Tracks))
deksummaries_merged_QC <- select(deksummaries_merged_QC, -c(
  `std.error_(Intercept)`, `std.error_Metadata_Time`, `statistic_(Intercept)`,
  `statistic_Metadata_Time`, `p.value_(Intercept)`, `p.value_Metadata_Time`))

# creates unique identifier for each well of each replicate
deksummaries_merged_QC$PlatesiRNAWell <- with(deksummaries_merged_QC,
    paste0( Metadata_Plate, "--", Metadata_siRNA, "--W", as.character(Metadata_WellNr)))
  
# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
# 3. Quality control on position summaries of all replicates
# adds quality control (QC) columns
# HTM_qc_intensity marks low nuclear intensity (here < 37.41) as FALSE
# HTM_qc_min marks positions with too few detected nuclei (here < 37.41) as FALSE
# HTM_qc_focus marks out of focus images (here 95QuantPowerLogLog -2.21 <= x <= 1.72) as FALSE
# HTM_qc_manual marks potential hits which escaped QC on a small margin but are 
# FALSE due to visual inspection
deksummaries_merged_QC %<>% mutate(HTM_qc_intensity = ifelse(
  mean_mean_integInt_nuclei < 37.41, FALSE, TRUE))
deksummaries_merged_QC %<>% mutate(HTM_qc_min = ifelse(min_nuclei < 48, FALSE, TRUE))
deksummaries_merged_QC %<>% mutate(HTM_qc_focus = ifelse(
  `95quant_PowerLogLog` < -2.21 | `95quant_PowerLogLog` >= -1.72, FALSE, TRUE))
deksummaries_merged_QC %<>% mutate(HTM_qc_death = ifelse(
  estimate_Metadata_Time < -0.5, FALSE, TRUE))
manual <- c("s19552", "s13490", "s4805", "s12362", "s8591", "s2838", "s30853", "s11727", "s2839",
  "s11719", "s11723", "s11736", "s11902", "s11904", "s11941", "s11949", "s14561", "s16491",
  "s17426", "s22035", "s22036", "s22120", "s22123", "s224035", "s32295", "s34942", "s194488",
  "s40362", "s6994")
deksummaries_merged_QC %<>% mutate(HTM_qc_manual = ifelse(
  Metadata_siRNA %in% manual, FALSE, TRUE))
# HTM_qc_cont marks wells with visible contamination
cont1 <- "May-Vogel-DEK-Bodies-Plate-01-batch1-01--s22356--W57"
cont2 <- "May-Vogel-DEK-Bodies-Plate-01-batch1-02--s223227--W59"
cont3 <- "May-Vogel-DEK-Bodies-Plate-01-batch1-02--s15830--W69"
cont4 <- "May-Vogel-DEK-Bodies-Plate-01-batch1-02--s2839--W82"
cont5 <- "May-Vogel-DEK-Bodies-Plate-03-batch1-01--s5020--W12"
cont6 <- "May-Vogel-DEK-Bodies-Plate-08-batch1-01--s14309--W31"
cont7 <- "May-Vogel-DEK-Bodies-Plate-06-batch1-02--s813--W23"
cont8 <- "May-Vogel-DEK-Bodies-Plate-04-batch1-02--s224035--W80"
cont_pos <- c(cont1, cont2, cont3, cont4, cont5, cont6, cont7, cont8)
deksummaries_merged_QC %<>% mutate(HTM_qc_cont = ifelse(
  PlatesiRNAWell %in% cont_pos, FALSE, TRUE))

# merges them together to a final column HTM_qc
deksummaries_merged_QC %<>% mutate(HTM_qc = ifelse(HTM_qc_min == FALSE |
  HTM_qc_intensity == FALSE | HTM_qc_focus == FALSE | HTM_qc_cont == FALSE | 
  HTM_qc_death == FALSE, FALSE, TRUE))

# when more than two positions within one well fail QC, 
# sets fourth position to fail as well, because data based on only one position is not robust
deksummaries_merged_QC %<>% group_by(PlatesiRNAWell) %>% mutate(HTM_qc =
  ifelse(sum(HTM_qc == TRUE) == 1 & HTM_qc == TRUE, FALSE, HTM_qc)) %>% ungroup()
# creates N_pos_siRNA column containing number of positions that passed QC per siRNA
deksummaries_merged_QC %<>% group_by(Metadata_siRNA) %>% mutate(N_pos_siRNA =
  sum(HTM_qc)) %>% ungroup()
# when less than 5 positions for the same siRNA pass QC, sets all positions for this siRNA to fail
deksummaries_merged_QC %<>% mutate(HTM_qc = ifelse(N_pos_siRNA < 5, FALSE, HTM_qc))
# creates N_pos column containing number of positions that passed QC per gene
deksummaries_merged_QC %<>% group_by(Metadata_treatment) %>% mutate(N_pos =
  sum(HTM_qc)) %>% ungroup()

# creates QC_type column describing the reason for QC fail
deksummaries_merged_QC$QC_type <- ""
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(
  HTM_qc_intensity == FALSE, paste(QC_type, "low DEK", sep = ";"), QC_type))
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(
  HTM_qc_min == FALSE, paste(QC_type, "low cells", sep = ";"), QC_type))
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(
  HTM_qc_death == FALSE, paste(QC_type, "cell death", sep = ";"), QC_type))
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(
  HTM_qc_focus == FALSE, paste(QC_type, "focus", sep = ";"), QC_type))
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(
  HTM_qc_cont == FALSE, paste(QC_type, "cont.", sep = ";"), QC_type))
deksummaries_merged_QC$QC_type <- sub(';', "", deksummaries_merged_QC$QC_type, perl = TRUE)
deksummaries_merged_QC %<>% mutate(QC_type = ifelse(QC_type == "", NA, QC_type))

#adds also a layout column
deksummaries_merged_QC$layout <- gsub("-Replicate\\d", "", deksummaries_merged_QC$PlateName)

# in case image paths have to be changed to local paths
deksummaries_merged_QC <- mutate(deksummaries_merged_QC,
  FullPath_Outlines = paste0( "D:/local_path/", Metadata_Plate, "/", gsub(".*\\/", "", FullPath_Outlines)),
  FullPath_Tracks = paste0("D:/local_path/", Metadata_Plate, "/", gsub(".*\\/", "", FullPath_Tracks)))

# save quality controlled position summaries
write_delim(deksummaries_merged_QC, paste0("tables/deksummaries_merged_QC.csv"), delim = ",")

## nucleisummaries ####
# concatenates several nucleisummaries.csv (must be in same directory) into a single dataframe
files <- list.files(path = "tables/", pattern = "nucleisummaries.csv", full.names = TRUE)
temp <- lapply(files, fread, sep = ",")
nucleisummaries_merged <- rbindlist(temp, fill = TRUE)

# save nuclei summaries
write_delim(nucleisummaries_merged, paste0("tables/nucleisummaries_merged.csv",), delim = ";")

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
# 4. Calculation of Z-scores to score hits

# log2 transformation of position summaries (deksummaries_merged_QC), 
# plate-wise calculation of negative control mean and subsequent subtraction of control from each 
# measurement

# filters deksummaries_merged_QC to contain only measurements of interest and metadata
measurements = c("median_tracklength_trackswdek", "mean_tracklength_trackswdek",
  "mean_ratio_onedek", "mean_ratio_twodek", "mean_ratio_threedek", "mean_ratio_fourdek", 
  "mean_ratio_fivedek", "mean_ratio_trackonedek", "mean_ratio_tracktwodek", 
  "mean_ratio_trackthreedek", "mean_ratio_trackfourdek", "mean_ratio_trackfivedek",
  "fraction_of_trackswmorethanonedek", "fraction_of_trackswmorethantwodek", 
  "fraction_of_trackswmorethanthreedek","fraction_of_trackswmorethanfourdek", 
  "fraction_of_trackswmorethanfivedek", "fraction_of_frameswmorethanonedek", 
  "fraction_of_frameswmorethantwodek", "fraction_of_frameswmorethanthreedek", 
  "fraction_of_frameswmorethanfourdek", "fraction_of_frameswmorethanfivedek", 
  "median_tracklifetime", "mean_tracklifetime", "mean_mean_dekbody_count", 
  "mean_mean_dekbody_count_seq", "mean_mean_intensity_dek", "fraction_of_frames_with_dek", 
  "fraction_of_tracks_with_dek", "mean_bodyduration", "mean_bodyduration_tracks", 
  "mean_mean_meanInt_nuclei", "mean_mean_integInt_nuclei", "mean_max_dek", "mean_max_dek_seq")
metadata = c("Metadata_siRNA", "PlateName", "Metadata_treatment", "Metadata_WellNr", 
  "Metadata_Pos", "replicate", "N_pos", "N_pos_siRNA", "RowCol", "HTM_qc", "layout")
  
# log2 transforms position measurements only for positions that passed quality control 
# to calculate means and sd of NEG1
log2_filt <- deksummaries_merged_QC %>% filter(HTM_qc == TRUE) %>% 
  select_(.dots = measurements) %>% log2()
# cleans up "-Inf" values from log2(0) 
is.na(log2_filt) <- log2_filt == "-Inf"
# and brings back the metadata
deksummaries_merged_filt_metadata <-
  deksummaries_merged_QC %>% filter(HTM_qc == TRUE) %>% select_(.dots = metadata)
log2_filt <- cbind(deksummaries_merged_filt_metadata, log2_filt)

# groups by gene name and plate to summarize negative control sNEG1 plate-wise, 
# calculates means and sd
mean_controls <- log2_filt %>% select(-c(Metadata_WellNr, Metadata_siRNA, Metadata_Pos,
  N_pos, RowCol, HTM_qc, layout)) %>% filter(Metadata_treatment == "sNEG1") %>%
  group_by(PlateName, replicate, Metadata_treatment) %>% 
  summarise_all(funs(mean = mean(., na.rm = TRUE))) %>% ungroup()
is.na(mean_controls) <- mean_controls == "NaN"
mean_controls %<>% select(-c(Metadata_treatment, replicate))
sd_controls <- log2_filt %>% select(-c(Metadata_WellNr, Metadata_siRNA, Metadata_Pos, N_pos, 
  RowCol, HTM_qc, layout)) %>% filter(Metadata_treatment=="sNEG1") %>% 
  group_by(PlateName, replicate, Metadata_treatment) %>% 
  summarise_all(funs(sd=sd(., na.rm=T))) %>% ungroup()
is.na(sd_controls) <- sd_controls == "NaN"
sd_controls %<>% select(-c(Metadata_treatment, replicate))

# log2 transforms position measurements for all positions
log2 <- deksummaries_merged_QC %>% select_(.dots = measurements) %>% log2()
is.na(log2) <- log2 == "-Inf"
deksummaries_merged_metadata <- deksummaries_merged_QC %>% select_(.dots = metadata)
log2 <- cbind(deksummaries_merged_metadata, log2)

# subtracts plate means of neg. control from each position measurement
metadata_siRNA <- unique(deksummaries_merged_metadata[, c('Metadata_siRNA', 'Metadata_treatment',
  'Metadata_WellNr', "layout", "N_pos", "N_pos_siRNA", "RowCol")]) %>% ungroup()
metadata_wowell <- unique(deksummaries_merged_metadata[, c('Metadata_siRNA', 'Metadata_treatment', 
  "PlateName", "replicate", "N_pos", "N_pos_siRNA")]) %>% ungroup()
metadata_worep <- unique(deksummaries_merged_metadata[, c('Metadata_siRNA', 'Metadata_treatment',
  'Metadata_WellNr', "layout", "N_pos", "N_pos_siRNA", "RowCol", "HTM_qc")]) %>% ungroup()
deksummaries_merged_mean_controls <- merge(deksummaries_merged_metadata, mean_controls)
deksummaries_merged_mean_controls <- select(deksummaries_merged_mean_controls,
    median_tracklength_trackswdek_mean:mean_max_dek_seq_mean)
log2_measurements <- select(log2, median_tracklength_trackswdek:mean_max_dek_seq)
substr_temp <- log2_measurements - deksummaries_merged_mean_controls

# divides each measurement by plate-wise SDs of neg control to calculate z-score
deksummaries_merged_sd_controls <- merge(deksummaries_merged_metadata, sd_controls)
deksummaries_merged_sd_controls <- select(deksummaries_merged_sd_controls,
    median_tracklength_trackswdek_sd:mean_max_dek_seq_sd)
zscore_temp <- substr_temp / deksummaries_merged_sd_controls
deksummaries_merged_zscore <- cbind(deksummaries_merged_metadata, zscore_temp)

# merges z-scores of single positions with unnormalized measurements from deksummaries_merged_QC
deksummaries <- deksummaries_merged_zscore %>% select(c(Metadata_siRNA, PlateName, 
  Metadata_WellNr, Metadata_Pos), measurements) %>% unique()
colnames(deksummaries) <- paste0(colnames(deksummaries), "_z_score")
deksummaries <- rename(deksummaries, Metadata_siRNA = Metadata_siRNA_z_score,
    PlateName = PlateName_z_score, Metadata_WellNr = Metadata_WellNr_z_score,
    Metadata_Pos = Metadata_Pos_z_score)
deksummaries <- merge(deksummaries, deksummaries_merged_QC, by = c(
  "Metadata_siRNA", "PlateName", "Metadata_WellNr", "Metadata_Pos"))

#creates siRNA summaries from unnormalized position summaries
siRNA_summaries_raw <- deksummaries_merged_QC %>% filter(HTM_qc == TRUE) %>% 
  select("Metadata_siRNA", measurements) %>% group_by(Metadata_siRNA) %>% 
  summarise_all(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE))) %>% ungroup()
is.na(siRNA_summaries_raw) <- siRNA_summaries_raw == "NaN"
colnames(siRNA_summaries_raw) <- paste0(colnames(siRNA_summaries_raw), "__siRNAraw")
siRNA_summaries_raw <- rename(siRNA_summaries_raw, Metadata_siRNA = Metadata_siRNA__siRNAraw)

# calculates different summaries from z-scores, regarding only positions that passed quality control
# groups plate-wise normalized measurements by siRNA to summarize each siRNA over all plates,
# calculates means and sd
siRNA_summaries_z <- deksummaries_merged_zscore %>% filter(HTM_qc == TRUE) %>% 
  select(-c(N_pos_siRNA, Metadata_treatment, PlateName, N_pos, replicate, Metadata_WellNr,
  RowCol, HTM_qc, Metadata_Pos, layout)) %>% group_by(Metadata_siRNA) %>% 
  summarise_all(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE))) %>% ungroup()
is.na(siRNA_summaries_z) <- siRNA_summaries_z == "NaN"
colnames(siRNA_summaries_z) <- paste0(colnames(siRNA_summaries_z), "__siRNA")
siRNA_summaries_z <- rename(siRNA_summaries_z, Metadata_siRNA = Metadata_siRNA__siRNA)

# groups plate-wise normalized measurements by gene to summarize each gene over all plates, 
# calculate means and sd
gene_summaries_z <- deksummaries_merged_zscore %>% filter(HTM_qc == TRUE) %>% 
  select(-c(Metadata_siRNA, PlateName, N_pos, replicate, Metadata_WellNr, RowCol, HTM_qc,
  Metadata_Pos, layout)) %>% group_by(Metadata_treatment) %>% 
  summarise_all(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE))) %>% ungroup()
is.na(gene_summaries_z) <- gene_summaries_z == "NaN"

# groups plate-wise normalized measurements by gene and replicate to summarize 
# each gene over each plate, calculates means and sd
treatmentsummaries_z <- deksummaries_merged_zscore %>% filter(HTM_qc == TRUE) %>% 
  select(-c(Metadata_siRNA, PlateName, N_pos, Metadata_WellNr, RowCol, HTM_qc, Metadata_Pos, 
  layout)) %>% group_by(Metadata_treatment, replicate) %>% 
  summarise_all(funs(mean = mean(., na.rm = TRUE), sd = sd(., na.rm = TRUE))) %>% ungroup()
is.na(treatmentsummaries_z) <- treatmentsummaries_z == "NaN"
colnames(treatmentsummaries_z) <- paste0(colnames(treatmentsummaries_z), "__rep")
treatmentsummaries_z <- rename(treatmentsummaries_z, 
  Metadata_treatment = Metadata_treatment__rep, replicate = replicate__rep)

# spreads siRNA summaries to allow for simultaneous plotting of several siRNAs per gene
siRNAmeans <- siRNA_summaries_z %>% select(ends_with("__siRNA")) %>% colnames()
treatmentsummaries_siRNAs <- merge(metadata_wowell, siRNA_summaries_z) %>% 
  select(-c(PlateName, N_pos, replicate)) %>% unique() %>% arrange(Metadata_treatment)
write_delim(treatmentsummaries_siRNAs, "tables/treatmentsummaries_siRNAs.csv", delim = ";")
# in excel: manually add column "siRNA" with identifier for unique siRNAs (1-3)
siRNA <- setDT(read_delim("tables/treatmentsummaries_siRNAs_excel.csv", delim = ";"))
siRNA <- select(siRNA, c(Metadata_treatment, Metadata_siRNA, siRNA))
treatmentsummaries_siRNAs <- merge(treatmentsummaries_siRNAs, siRNA)
treatmentsummaries_siRNAs <- merge(treatmentsummaries_siRNAs, siRNA_summaries_z) %>% setDT()
treatmentsummaries_siRNAs <- melt(treatmentsummaries_siRNAs, measure.vars = patterns("__siRNA"))
treatmentsummaries_siRNAs <- dcast(treatmentsummaries_siRNAs,
  Metadata_treatment ~ variable + siRNA, fun = mean, value.vars = "value")
is.na(treatmentsummaries_siRNAs) <- treatmentsummaries_siRNAs == "NaN"

# spreads treatment summaries to allow for simultaneous plotting of plate replicate values of genes
replicatemeans <- treatmentsummaries_z %>% select(ends_with("__rep")) %>% colnames()
setDT(treatmentsummaries_z)
treatmentsummaries_treatments <- melt(treatmentsummaries_z, measure.vars = patterns("__rep"))
treatmentsummaries_treatments <- dcast(treatmentsummaries_treatments,
  Metadata_treatment ~ variable + replicate, fun = mean, value.vars = "value")
is.na(treatmentsummaries_treatments) <- treatmentsummaries_treatments == "NaN"

# merge all z-score summaries together
treatmentsummaries <- merge(metadata_worep, gene_summaries_z)
treatmentsummaries <- merge(treatmentsummaries, treatmentsummaries_siRNAs)
treatmentsummaries <- merge(treatmentsummaries, treatmentsummaries_treatments)
siRNA_summaries <- merge(metadata_siRNA, siRNA_summaries_z) 
siRNA_summaries <- merge(siRNA_summaries, siRNA) %>% 
  select(Metadata_siRNA, Metadata_treatment, siRNA, everything())

# creates data frame with selected measurements and tidies up for excel
excel <- treatmentsummaries %>% select(matches(
      'ratio_trackonedek | ratio_tracktwodek | ratio_trackthreedek | ratiotrackfourdek |
      intensity_dek|dekbody_count|max_dek'))
excel_metadata <- treatmentsummaries %>% select_(.dots = names(metadata_worep))
excel <- cbind(excel_metadata, excel)
excel_comb <- excel %>% filter(HTM_qc == TRUE) %>% group_by(Metadata_treatment) %>% 
  summarise(layout = paste(layout, collapse = ","), 
  Metadata_WellNr = paste(Metadata_WellNr, collapse=","), 
  RowCol = paste(RowCol, collapse = ","), 
  Metadata_siRNA = paste(Metadata_siRNA, collapse = ",")) %>% ungroup()
excel %<>% select(-c(HTM_qc, Metadata_siRNA, Metadata_WellNr, RowCol, layout)) %>% unique()
excel <- merge(excel_comb, excel) %>% rename(gene = Metadata_treatment, siRNA = Metadata_siRNA,
  well = Metadata_WellNr)
excel %<>% select(c(gene, siRNA, layout, well, RowCol, N_pos), order(colnames(.))) %>% 
  select(-N_pos_siRNA) %>% unique()

# creates a gene list of genes that passed QC to use as input for GO enrichment
input <- treatmentsummaries %>% filter(HTM_qc == TRUE & N_pos_siRNA > 6 & 
  !Metadata_treatment %in% c("empty", "sNEG1")) %>% select(Metadata_treatment) %>% 
  unique() %>% t() %>% c()

#saves all z-score summaries
write.table(input, "GOinput.txt", row.names = F, quote = FALSE)
write_delim(deksummaries, "tables/deksummaries.csv", delim = ";")
write_delim(treatmentsummaries, "tables/treatmentsummaries.csv", delim =";")
write_delim(treatmentsummaries_siRNAs, "tables/treatmentsummaries_siRNAs.csv", delim = ";")
write_delim(excel, "tables/treatmentsummaries_excel.csv", delim = ";")
write_delim(siRNA_summaries, "tables/siRNAsummaries.csv", delim = ";")







