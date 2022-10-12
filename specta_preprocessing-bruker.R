library("MALDIquant")
library("MALDIquantForeign")

#Three arguments are required: input directory of .mzml files == args[1], and output directory were .csv files of peaks should go == args[2] 
args = commandArgs(trailingOnly=TRUE)

#import unprocessed mzml files, exported from Launchpad
spectra <- import(args[1], removeEmptySpectra = TRUE)

spectra_processed <- transformIntensity(spectra,method="sqrt")
spectra_processed  <- trim(spectra_processed, range = c(4000,19090))
spectra_processed  <- smoothIntensity(spectra_processed , method="SavitzkyGolay",halfWindowSize=20)
rm(spectra)
spectra_bl_removed  <- removeBaseline(spectra_processed, method="SNIP", iterations=40)
rm(spectra_processed)
spectra_cal <- calibrateIntensity(spectra_bl_removed, method="median")

peaks2 <- detectPeaks(spectra_cal, method="SuperSmoother", halfWindowSize= 20, SNR=2)


# to calibrate
ref_masses = c(4364.4, 5095.7619, 6371.5031, 6446.3097, 6541.7186,7273.3643, 7288.8476, 8499.85469999999, 9006.383, 9704.33379999999, 10430.1588, 11564.2018999999, 11580.351, 11735.440999, 12769.46, 13133.0696, 13540.873,14126.3943999999, 14875.2147999999, 15281.0238, 15768.8566999999,17603.1810999999, 17711.3859999999)
ref_massesH = ref_masses +1
refPeaks <- createMassPeaks(mass=ref_massesH, intensity=rep(1, length(ref_massesH)))

#warp spectra
warpedPeaks <- list()
warpedSpectra <- list()
for (i in (1:length(peaks2))){
  warpingFunctions <- try(determineWarpingFunctions(peaks2[i], reference=refPeaks, method="quadratic",tolerance = 0.002, plot=FALSE, plotInteractive=TRUE))
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
