library("MALDIquant")
library("MALDIquantForeign")

#Three arguments are required: input directory of .mzml files == args[1], and output directory were .csv files of peaks should go == args[2] and a csv file which contains the positions on the MALDI targes as well as the TGNR == args[3] ()
args = commandArgs(trailingOnly=TRUE)

spectra_processed <- transformIntensity(spectra,method="sqrt")
spectra_processed  <- trim(spectra_processed, range = c(4000,19090))
spectra_processed  <- smoothIntensity(spectra_processed , method="SavitzkyGolay",halfWindowSize=20)
rm(spectra)
spectra_bl_removed  <- removeBaseline(spectra_processed, method="SNIP", iterations=160)
rm(spectra_processed)
spectra_cal <- calibrateIntensity(spectra_bl_removed, method="median")

peaks2 <- detectPeaks(spectra_cal, method="SuperSmoother", halfWindowSize= 20, SNR=2)
plot(spectra_cal[[1]])
library(gridExtra)
library(grid)
library(ggplot2)


## label highest peaks (top 10) and avoid label overlap
plot(spectra_cal[[1]])
top10 <- intensity(peaks2[[1]]) %in% sort(intensity(peaks2[[1]]), decreasing=TRUE)[1:10]
labelPeaks(peaks2[[1]], index=top10)


# to calibrate
ref_masses = c(4364.4, 5095.7619, 6371.5031, 6446.3097, 6541.7186,7273.3643, 7288.8476, 8499.85469999999, 9006.383, 9704.33379999999, 10430.1588, 11564.2018999999, 11580.351, 11735.440999, 12769.46, 13133.0696, 13540.873,14126.3943999999, 14875.2147999999, 15281.0238, 15768.8566999999,17603.1810999999, 17711.3859999999)
ref_massesH = ref_masses +1
refPeaks <- createMassPeaks(mass=ref_massesH, intensity=rep(1, length(ref_massesH)))

#warp spectra
warpedPeaks <- list()
warpedSpectra <- list()
for (i in (1:length(peaks2))){
  warpingFunctions <- try(determineWarpingFunctions(peaks2[i], reference=refPeaks, method="linear",tolerance = 0.002, plot=FALSE, plotInteractive=TRUE))
  warpedPeaks[i] <-try(warpMassPeaks(peaks2[i], warpingFunctions))
  warpedSpectra[i] <- try(warpMassSpectra(spectra_cal[i], warpingFunctions))
}

#Find indices of files for which the warping did not work
warping_worked<-list()
for (i in (1:length(warpedSpectra))){
  warping_worked[i]<-class(warpedSpectra[[i]])!='character'
}

warping_worked<-unlist(warping_worked)

#include only files for which warpping did work
warpedSpectra <- warpedSpectra[warping_worked]
warpedPeaks<- warpedPeaks[warping_worked]

#extract files
files<-list()
for (i in (1:length(warpedPeaks))){
  files[i] <- list(sapply(warpedPeaks[i],function(x)metaData(x)$file))
}

files<-as.data.frame(unlist(files))
colnames(files)<-'files'

# Add positions to separate column to merge and remove '.mzml'
files$rppos<-gsub('.*\\/','',files[,'files'])
files$rppos<-gsub('.mzXML', '', files$rppos)

#read position
poso.tgnr.invas<-read.csv2(args[3], sep=',')
#Harmonize naming of run
poso.tgnr.invas$rppos<-gsub('mabr_D19_0510_','mabr_d19_0510_', poso.tgnr.invas$rppos)

#merge 'files' and poso.tgnr.invas
files.tgnr<-merge(files, poso.tgnr.invas, by='rppos', all.x = TRUE, all.y= FALSE)
#files.tgnr['new_files']<-paste0(gsub('[^\\/]*mzXML$','',files.tgnr$files),files.tgnr$TGNR,'.mzXML')
files.tgnr['new_files']<-paste0(gsub('[^\\/]*mzXML$','',files.tgnr$files),files.tgnr$TGNR,'_', gsub('.+\\_','',files.tgnr$rppos),'.mzXML')


#rename file in Peaks object, that it is exported as csv with the TGNR as name
for (i in (1:length(warpedPeaks))){
  for (j in (1:length(files.tgnr$files))){
       warpedPeaks[[i]]@metaData$file<-ifelse(warpedPeaks[[i]]@metaData$file == files.tgnr$files[j] & !is.na(files.tgnr$TGNR[j]), files.tgnr$new_files[j],warpedPeaks[[i]]@metaData$file)
  }
}


export(warpedPeaks, type="csv", path=args[2], force=TRUE)

#summarise all peaks into dataframe
peaks_df<-data.frame()
for (i in (1:length(warpedPeaks))){
  peaks_df_temp<-data.frame(warpedPeaks[[i]]@mass)
  peaks_df_temp<-cbind(peaks_df_temp, warpedPeaks[[i]]@intensity)
  peaks_df_temp<-cbind(peaks_df_temp, rep(gsub('(^.+\\/)(.*)(\\.mzXML$)','\\2', warpedPeaks[[i]]@metaData$file), length(warpedPeaks[[i]])))
  peaks_df<-rbind(peaks_df, peaks_df_temp)
}

write.csv(peaks_df, paste0(args[2], 'peaks_df.csv'))
