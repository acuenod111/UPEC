# plot peak occurences by phylogroup
# library
library('MALDIquant')
library('MALDIquantForeign')
library('tidyr')
library('cowplot')
library('ggplot2')
library('dplyr')

# Iimport MALDIQuant picked peaks
shimadzu <- read.csv2('/Users/aline/Doc.Mobility/01_Data/01_spectra/Shimadzu_csv_medianpeaks_df.csv', sep=',')
colnames(shimadzu) <- c('mass', 'intensity', 'spectrum')
shimadzu['TGNR']<-gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', shimadzu$spectrum)

# Bruker
Bruker <- read.csv2('/Users/aline/Doc.Mobility/01_Data/01_spectra/Bruker_csv_medianpeaks_df.csv', sep = ',')
colnames(Bruker) <- c('mass', 'intensity', 'spectrum')
Bruker['code']<-gsub('([[:alnum:]]{8}\\-[[:alnum:]]{4}\\-[[:alnum:]]{4}\\-[[[:alnum:]]{4}\\-[[[:alnum:]]{12})(\\..*$)', '\\1', Bruker$spectrum)


# translate where necessary
trans <- read.csv2('/Users/aline/Doc.Mobility/01_Data/01_spectra/E.coli/Bruker/Brukercode.tgnr.csv', header = T, sep = ',')
colnames(trans)<-c('code', 'TGNR')
Bruker <- merge(Bruker, trans, by = 'code', all = T, all.y = F)
Bruker['TGNR']<-ifelse(is.na(Bruker$TGNR), gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', Bruker$spectrum), Bruker$TGNR)
Bruker$code <- NULL

# import dereplicated strains (one strain per clinical case (n=828))
derep<-read.csv2('/Users/aline/Doc.Mobility/01_Data/05_Strains/one_per_case_strain_828.txt', header = F)
colnames(derep)<-'fullpath'
derep['assembly'] <- gsub('(.*\\/)(.*)', '\\2', derep$fullpath)
derep['strain']<-gsub('\\.fna', '', derep$assembly)
derep['TGNR']<-gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', derep$strain)

# only consider dereplicated strains
# subset shimadzu
shimadzu <- shimadzu[shimadzu$TGNR %in% derep$TGNR,]
shimadzu <- shimadzu[shimadzu$TGNR %in% Bruker$TGNR,]
shimadzu$mass<-as.numeric(shimadzu$mass)
shimadzu$intensity<-as.numeric(shimadzu$intensity)

# subset Bruker files
Bruker <- Bruker[Bruker$TGNR %in% derep$TGNR,]
Bruker <- Bruker[Bruker$TGNR %in% shimadzu$TGNR,]
Bruker$mass<-as.numeric(Bruker$mass)
Bruker$intensity<-as.numeric(Bruker$intensity)

# build mass peaks object in order to build intensity matrix
# shimadzu
shimadzu_peaks <- list()
for (spectrum in unique(shimadzu$spectrum)){
  shimadzu_peaks[[spectrum]]<-createMassPeaks(shimadzu[shimadzu$spectrum == spectrum, 'mass'], shimadzu[shimadzu$spectrum == spectrum, 'intensity'])
}
shimadzu_peaks <- binPeaks(shimadzu_peaks, tolerance = 500 / 1000000, method = "relaxed")
shimadzu_peaks <- filterPeaks(shimadzu_peaks, minFrequency=0.10)
shimadzu_intensity_matrix <- intensityMatrix(shimadzu_peaks)
# set all non-present peaks to intensity 0
shimadzu_intensity_matrix[is.na(shimadzu_intensity_matrix)] <- 0

rownames(shimadzu_intensity_matrix)<-names(shimadzu_peaks)

# visually inspect binary peaks
featureMatrix_binary <- shimadzu_intensity_matrix
featureMatrix_binary[!featureMatrix_binary == 0] <- 1
# remove all which are always or never there
featureMatrix_binary <- as.data.frame(featureMatrix_binary)
featureMatrix_binary <- featureMatrix_binary %>% select(where(~length(unique(.)) > 1))
featureMatrix_binary['spectrum']<-rownames(featureMatrix_binary)
featureMatrix_binary['TGNR']<-gsub('\\_.*$', '', featureMatrix_binary$spectrum)

# import phylogroup
phylo<-read.delim('/Users/aline/Doc.Mobility/01_Data/03_sequences/14_phylogroups/mash_phylo.csv',quote="", sep=';',header=T, fill = TRUE)
colnames(phylo)<-c("sample","best_mash_dist","Phylogroup","medoid.phylogroup")
# only look at dereplicated strains
phylo['label'] <- phylo$sample
phylo <- phylo[phylo$label %in% derep$strain,]
phylo['TGNR']<-gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', phylo$label)

featureMatrix_binary<-merge(featureMatrix_binary, phylo[,c('TGNR', 'Phylogroup')], by = 'TGNR', all.x = T, all.y = F)

table(featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR), 'Phylogroup'])
# remove unfrequent phylogroups
featureMatrix_binary <- featureMatrix_binary[!featureMatrix_binary$Phylogroup %in% c('cladeV', 'E', 'E1','E2','G','Unknown'),]

# summarise by phylogroup
featureMatrix_binary_sum <- featureMatrix_binary %>% 
  group_by(Phylogroup) %>% 
  summarise(across(colnames(featureMatrix_binary)[!colnames(featureMatrix_binary) %in% c('TGNR', 'Phylogroup', 'spectrum')], mean))
table(featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR), 'Phylogroup'])

# fin peaks which are not frequent (less than 25%) in one cluster
binary_sum_sel <- featureMatrix_binary_sum[,apply(featureMatrix_binary_sum, 2, function(X) any(X <= 0.25|grepl('[[:alpha:]]', X))), drop=FALSE]
# find peaks which are present in at least 75% of the other group
binary_sum_sel <- binary_sum_sel[!is.na(featureMatrix_binary_sum$Phylogroup),apply(binary_sum_sel, 2, function(X) any(X > 0.5|grepl('[[:alpha:]]', X))), drop=FALSE]

# In order to draw a heatmap of occurences, convert to 'long' format
# featureMatrix_binary_sum_long<-pivot_longer(featureMatrix_binary_sum, cols = colnames(featureMatrix_binary_sum)[!colnames(featureMatrix_binary_sum) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')
featureMatrix_binary_sum_long<-pivot_longer(binary_sum_sel, cols = colnames(binary_sum_sel)[!colnames(binary_sum_sel) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')

# order according to increasing masses
featureMatrix_binary_sum_long$mass<-factor(featureMatrix_binary_sum_long$mass, levels = as.character(sort(unique(as.numeric(featureMatrix_binary_sum_long$mass)))))
# remove unfrequent phylogroups
featureMatrix_binary_sum_long <- featureMatrix_binary_sum_long[!featureMatrix_binary_sum_long$Phylogroup %in% c(NA, "E1","E2", "G"),]
featureMatrix_binary_sum_long$Phylogroup <- factor(featureMatrix_binary_sum_long$Phylogroup, levels = c("F", "B2-1", "B2-2","D1","D2","D3","A", "C", "B1"))
featureMatrix_binary_sum_long_shimadzu <- featureMatrix_binary_sum_long


# plot all masses
occurence_plot_shimadzu_phylo <- ggplot(featureMatrix_binary_sum_long, aes( x= mass, y = Phylogroup, fill=occurence)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('')


# import PAPG VARIANT
pap<-read.delim('/Users/aline/Doc.Mobility/01_Data/03_sequences/06_nonribosomal-targets/pagG_var.csv',quote="", sep=';',header=T, fill = TRUE)
pap$X.FILE <- gsub("103713-19.fna", "103713-49.fna", pap$X.FILE)
pap <- pap[pap$X.FILE %in% derep$assembly,]
pap['TGNR']<-gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', pap$strain)
# summarise all papGII containing ones
pap$papG_variant <- gsub('no papG_', '',pap$papG_variant)
pap$papG_variant <- ifelse(grepl('papGII\\_', pap$papG_variant), 'papGII', pap$papG_variant)

# summarise all papG variants other than papGII to 'Other papG'
pap$papG_variant <- gsub("papGIII|papGIV|papGV", "Other papG", pap$papG_variant)

featureMatrix_binary<-merge(featureMatrix_binary, pap[,c('TGNR', 'papG_variant')], by = 'TGNR', all.x = T, all.y = F)
featureMatrix_binary_check<-featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR),]
table(featureMatrix_binary_check$Phylogroup, featureMatrix_binary_check$papG_variant)

# summarise by papG variant
featureMatrix_binary_sum <- featureMatrix_binary %>% 
  group_by(papG_variant) %>% 
  summarise(across(colnames(featureMatrix_binary)[!colnames(featureMatrix_binary) %in% c('TGNR', 'Phylogroup','papG_variant', 'spectrum')], mean))

# check how many per variant
table(featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR),'papG_variant'])

# fin peaks which are not frequent (less than 25%) in one cluster
binary_sum_sel <- featureMatrix_binary_sum[!is.na(featureMatrix_binary_sum$papG_variant),apply(featureMatrix_binary_sum, 2, function(X) any(X <= 0.25|grepl('[[:alpha:]]', X))), drop=FALSE]
# find peaks which are present in at least 75% of the other group
binary_sum_sel <- binary_sum_sel[,apply(binary_sum_sel, 2, function(X) any(X > 0.3|grepl('[[:alpha:]]', X))), drop=FALSE]


# In order to draw a heatmap of occurences, convert to 'long' format
# featureMatrix_binary_sum_long<-pivot_longer(featureMatrix_binary_sum, cols = colnames(featureMatrix_binary_sum)[!colnames(featureMatrix_binary_sum) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')
featureMatrix_binary_sum_long<-pivot_longer(binary_sum_sel, cols = colnames(binary_sum_sel)[!colnames(binary_sum_sel) %in% c('papG_variant','Phylogroup')],  names_to = "mass", values_to = 'occurence')

# order according to increasing masses
featureMatrix_binary_sum_long$mass<-factor(round(as.numeric(featureMatrix_binary_sum_long$mass)), levels = as.character(sort(unique(round(as.numeric(featureMatrix_binary_sum_long$mass))))))

# remove NA
featureMatrix_binary_sum_long <- featureMatrix_binary_sum_long[!is.na(featureMatrix_binary_sum_long$papG_variant),]

# relevel
featureMatrix_binary_sum_long$papG_variant <- factor(featureMatrix_binary_sum_long$papG_variant, levels = c("papGII", "Other papG","no papG"))
featureMatrix_binary_sum_long_shimadzu_pap <- featureMatrix_binary_sum_long

# plot all masses
occurence_plot_shimadzu_pap <- ggplot(featureMatrix_binary_sum_long_shimadzu_pap, aes(x= mass, y = papG_variant, fill=occurence)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('')

occurence_plot_shimadzu_pap


# Bruker
bruker_peaks <- list()
for (spectrum in unique(Bruker$spectrum)){
  bruker_peaks[[spectrum]]<-createMassPeaks(Bruker[Bruker$spectrum == spectrum, 'mass'], Bruker[Bruker$spectrum == spectrum, 'intensity'])
}
bruker_peaks <- binPeaks(bruker_peaks, tolerance = 500 / 1000000, method = "relaxed")
bruker_peaks <- filterPeaks(bruker_peaks, minFrequency=0.10)
bruker_intensity_matrix <- intensityMatrix(bruker_peaks)
bruker_intensity_matrix[is.na(bruker_intensity_matrix)] <- 0

rownames(bruker_intensity_matrix)<-names(bruker_peaks)

# visually inspect binary peaks
featureMatrix_binary <- bruker_intensity_matrix
featureMatrix_binary[!featureMatrix_binary == 0] <- 1
# remove all which are always or never there
featureMatrix_binary <- as.data.frame(featureMatrix_binary)
featureMatrix_binary <- featureMatrix_binary %>% select(where(~length(unique(.)) > 1))

featureMatrix_binary['spectrum']<-rownames(featureMatrix_binary)

featureMatrix_binary <- merge(featureMatrix_binary, Bruker[!duplicated(Bruker$spectrum),c("spectrum","TGNR")], by = 'spectrum')
featureMatrix_binary$TGNR <- gsub("120706-1-18-1", "120706-1-18", featureMatrix_binary$TGNR)
featureMatrix_binary<-merge(featureMatrix_binary, phylo[,c('TGNR', 'Phylogroup')], by = 'TGNR', all.x = T, all.y = F)
featureMatrix_binary_bruker <- featureMatrix_binary
# remove unfrequent phylogroups
featureMatrix_binary <- featureMatrix_binary[!featureMatrix_binary$Phylogroup %in% c('cladeV', 'E', 'G','Unknown'),]


featureMatrix_binary_sum <- featureMatrix_binary %>% 
  group_by(Phylogroup) %>% 
  summarise(across(colnames(featureMatrix_binary)[!colnames(featureMatrix_binary) %in% c('TGNR', 'Phylogroup', 'spectrum')], mean))
table(featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR), 'Phylogroup'])

# fin peaks which are not frequent (less than 25%) in one cluster
binary_sum_sel <- featureMatrix_binary_sum[!is.na(featureMatrix_binary_sum$Phylogroup),apply(featureMatrix_binary_sum, 2, function(X) any(X <= 0.25|grepl('[[:alpha:]]', X))), drop=FALSE]
# find peaks which are present in at least 75% of the other group
binary_sum_sel <- binary_sum_sel[,apply(binary_sum_sel, 2, function(X) any(X > 0.5|grepl('[[:alpha:]]', X))), drop=FALSE]

# In order to draw a heatmap of occurences, convert to 'long' format
# featureMatrix_binary_sum_long<-pivot_longer(featureMatrix_binary_sum, cols = colnames(featureMatrix_binary_sum)[!colnames(featureMatrix_binary_sum) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')
featureMatrix_binary_sum_long<-pivot_longer(binary_sum_sel, cols = colnames(binary_sum_sel)[!colnames(binary_sum_sel) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')

# order according to increasing masses
featureMatrix_binary_sum_long <- featureMatrix_binary_sum_long[!featureMatrix_binary_sum_long$Phylogroup %in% c(NA, "E1","E2", "G"),]
featureMatrix_binary_sum_long$Phylogroup <- factor(featureMatrix_binary_sum_long$Phylogroup, levels = c("F", "B2-1", "B2-2","D1","D2","D3","A", "C", "B1"))

featureMatrix_binary_sum_long_bruker <- featureMatrix_binary_sum_long

# round all masses and make levels, such that these are the same for shimadzu and bruker
featureMatrix_binary_sum_long_bruker$mass <- round(as.numeric(featureMatrix_binary_sum_long_bruker$mass))
featureMatrix_binary_sum_long_shimadzu$mass <- round(as.numeric(as.character(featureMatrix_binary_sum_long_shimadzu$mass)))
featureMatrix_binary_sum_long_bruker$mass<-factor(featureMatrix_binary_sum_long_bruker$mass, levels = as.character(sort(unique(as.numeric(featureMatrix_binary_sum_long_bruker$mass)))))

# shimadzu as well 
featureMatrix_binary_sum_long_shimadzu$mass<-factor(featureMatrix_binary_sum_long_shimadzu$mass, levels = as.character(sort(unique(as.numeric(featureMatrix_binary_sum_long_shimadzu$mass)))))

# plot shimadu
occurence_plot_shimadzu <- ggplot(featureMatrix_binary_sum_long_shimadzu, aes( x= mass, y = Phylogroup, fill=occurence)) + 
  geom_tile() +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('')


# plot bruker
occurence_plot_bruker <- ggplot(featureMatrix_binary_sum_long_bruker, aes( x= mass, y = Phylogroup, fill=occurence)) + 
  geom_tile() +
  scale_x_discrete(drop=FALSE) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('')

cowplot::plot_grid(occurence_plot_shimadzu + ggtitle('Shimadzu') +  theme(legend.position = 'none'), occurence_plot_bruker + ggtitle('Microflex') +  theme(legend.position = 'bottom'), nrow = 2, rel_heights = c(0.8,1))

pdf('/Users/aline/Doc.Mobility/04_Presentations/occurence_plot_peaks_bruker_phylogroup.pdf', width = 10, height = 4)
occurence_plot_bruker
dev.off()

pdf('/Users/aline/Doc.Mobility/04_Presentations/occurence_plot_peaks_shimadzu_phylogroup.pdf', width = 10, height = 4)
occurence_plot_shimadzu
dev.off()



pdf('/Users/aline/Doc.Mobility/04_Presentations/occurence_plot_peaks_phylogroup.pdf', width = 20, height = 8)
cowplot::plot_grid(occurence_plot_shimadzu + ggtitle('Shimadzu')+  theme(legend.position = 'bottom'), occurence_plot_bruker + ggtitle('Microflex') +  theme(legend.position = 'bottom'), nrow = 2, rel_heights = c(0.8,1))
dev.off()



# papG bruker
featureMatrix_binary<-merge(featureMatrix_binary_bruker, pap[,c('TGNR', 'papG_variant')], by = 'TGNR', all.x = T, all.y = F)
featureMatrix_binary_check<-featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR),]
table(featureMatrix_binary_check$Phylogroup, featureMatrix_binary_check$papG_variant)

# summarise by papG variant
featureMatrix_binary_sum <- featureMatrix_binary %>% 
  group_by(papG_variant) %>% 
  summarise(across(colnames(featureMatrix_binary)[!colnames(featureMatrix_binary) %in% c('TGNR', 'Phylogroup','papG_variant', 'spectrum')], mean))

# check how many per variant
table(featureMatrix_binary[!duplicated(featureMatrix_binary$TGNR),'papG_variant'])

# find peaks which are not frequent (less than 25%) in one cluster
binary_sum_sel <- featureMatrix_binary_sum[!is.na(featureMatrix_binary_sum$papG_variant),apply(featureMatrix_binary_sum, 2, function(X) any(X <= 0.25|grepl('[[:alpha:]]', X))), drop=FALSE]
# find peaks which are present in at least 75% of the other group
binary_sum_sel <- binary_sum_sel[!is.na(binary_sum_sel$papG_variant),apply(binary_sum_sel, 2, function(X) any(X > 0.3|grepl('[[:alpha:]]', X))), drop=FALSE]

# In order to draw a heatmap of occurences, convert to 'long' format
# featureMatrix_binary_sum_long<-pivot_longer(featureMatrix_binary_sum, cols = colnames(featureMatrix_binary_sum)[!colnames(featureMatrix_binary_sum) %in% 'Phylogroup'],  names_to = "mass", values_to = 'occurence')
featureMatrix_binary_sum_long<-pivot_longer(binary_sum_sel, cols = colnames(binary_sum_sel)[!colnames(binary_sum_sel) %in% c('papG_variant','Phylogroup')],  names_to = "mass", values_to = 'occurence')

# order according to increasing masses
featureMatrix_binary_sum_long$mass<-factor(round(as.numeric(featureMatrix_binary_sum_long$mass)), levels = as.character(sort(unique(round(as.numeric(featureMatrix_binary_sum_long$mass))))))

# remove NA
featureMatrix_binary_sum_long <- featureMatrix_binary_sum_long[!is.na(featureMatrix_binary_sum_long$papG_variant),]

# relevel papG variants
featureMatrix_binary_sum_long$papG_variant <- factor(featureMatrix_binary_sum_long$papG_variant, levels = c("papGII", "Other papG","no papG"))
featureMatrix_binary_sum_long_shimadzu <- featureMatrix_binary_sum_long

# plot all masses
occurence_plot_bruker_pap <- ggplot(featureMatrix_binary_sum_long, aes(x= mass, y = papG_variant, fill=occurence)) + 
  geom_tile() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient(low = 'white', high = 'black') + 
  xlab('') +
  ylab('')

occurence_plot_bruker_pap

cowplot::plot_grid(occurence_plot_bruker_pap, occurence_plot_shimadzu_pap+theme(legend.position = 'none'), ncol = 1)


pdf('/Users/aline/Doc.Mobility/04_Presentations/MALDI_pap.pdf', height=3.5, width =9)
cowplot::plot_grid(occurence_plot_bruker_pap+theme(legend.position = 'bottom'), occurence_plot_shimadzu_pap+theme(legend.position = 'bottom'), ncol = 2, rel_widths = c(1,2), labels = c('microflex Biotyper', 'Axima Confidence'))
dev.off()


