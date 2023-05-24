# visualise msa from PCR genes
# load packages
library(ggplot2)
library(ggstance)
library(dplyr)
library(tidyr)
library(stringr)
library(ggtree)
library(phytools)
library(phylotools)
library(bioseq)
library(seqinr)
library(ggmsa)
library(ggtreeExtra)
library(ggh4x)
library(Biostrings)
library(cowplot)

# gapC
# import gapC MSA
gapC_alignment <- readDNAStringSet("gapC.core.full.aln", format="fasta")
gapC_alignment <- DNAMultipleAlignment(gapC_alignment)

# locate primers
library(stringr)
gapC_ref <- as.character(gapC_alignment@unmasked$Reference)
#forward primer
str_locate_all(pattern ='CGCGGCAGAAAATATCATTCCC', gapC_ref)
primer_F_gapC <- tidy_msa(gapC_alignment, 597, 618)
#probe
str_locate_all(pattern ='ACGCGTGCCGGTGAAAACAGG', gapC_ref) #
probe_gapC <- tidy_msa(gapC_alignment, 693, 713)
#loop3
str_locate_all(pattern ='GGTCACTGAACTGGTATCGATTC', gapC_ref) # reverse complement of primer!
primer_R_gapC <- tidy_msa(gapC_alignment, 717, 739)

# summarise which variant of loops occur how frequently
primer_F_gapC_alignment_wide <- pivot_wider(primer_F_gapC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_F_gapC',2:(max(primer_F_gapC$position) - min(primer_F_gapC$position)+2), sep = '') 
probe_gapC_alignment_wide <- pivot_wider(probe_gapC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('probe_gapC',2:(max(probe_gapC$position) - min(probe_gapC$position)+2), sep = '') 
primer_R_gapC_alignment_wide <- pivot_wider(primer_R_gapC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_R_gapC',2:(max(primer_R_gapC$position) - min(primer_R_gapC$position)+2), sep = '') 

# merge both loops
gapC_alignments <- merge(primer_F_gapC_alignment_wide, probe_gapC_alignment_wide, by = 'name')
gapC_alignments <- merge(gapC_alignments, primer_R_gapC_alignment_wide, by = 'name')

# check which occur how often
#remove reference
gapC_alignments <- gapC_alignments[gapC_alignments$name != 'Reference',]
#check forward primer
count_primer_F_gapC <- primer_F_gapC_alignment_wide %>% 
  group_by(primer_F_gapC) %>%
  summarise(n = n())
count_primer_F_gapC[!grepl('N',count_primer_F_gapC$primer_F_gapC),] 
'CGCGGCAGAAAATATCATTCCC'=='CGCGGCAGAAAATATCATTCCC' #this corresponds to the primer used!

# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_F_gapC_alignment_wide_sum <- primer_F_gapC_alignment_wide
primer_F_gapC_alignment_wide_sum$primer_F_gapC <- ifelse(grepl('N',primer_F_gapC_alignment_wide_sum$primer_F_gapC), 'ambiguous sequence containing "N"', primer_F_gapC_alignment_wide_sum$primer_F_gapC)

gapC_F_plot <- ggplot(primer_F_gapC_alignment_wide_sum, aes(x = reorder(primer_F_gapC,primer_F_gapC,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('gapC forward primer')
gapC_F_plot

#check reverse primer
count_primer_R_gapC <- primer_R_gapC_alignment_wide %>% 
  group_by(primer_R_gapC) %>%
  summarise(n = n())
count_primer_R_gapC[!grepl('N',count_primer_R_gapC$primer_R_gapC),] 
'GGTCACTGAACTGGTATCGATTC'=='GGTCACTGAACTGGTATCGATTC' #this corresponds to the primer used!

# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_R_gapC_alignment_wide_sum <- primer_R_gapC_alignment_wide
primer_R_gapC_alignment_wide_sum$primer_R_gapC <- ifelse(grepl('N',primer_R_gapC_alignment_wide_sum$primer_R_gapC), 'ambiguous sequence containing "N"', primer_R_gapC_alignment_wide_sum$primer_R_gapC)

gapC_R_plot <- ggplot(primer_R_gapC_alignment_wide_sum, aes(x = reorder(primer_R_gapC,primer_R_gapC,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('gapC reverse primer')
gapC_R_plot

#check probe
count_probe_gapC <- probe_gapC_alignment_wide %>% 
  group_by(probe_gapC) %>%
  summarise(n = n())
count_probe_gapC[!grepl('N',count_probe_gapC$probe_gapC),] 
'ACGCGTGCCGGTGAAAACAGG'=='ACGCGTGCCGGTGAAAACAGG' #this corresponds to the probe used!


# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
probe_gapC_alignment_wide_sum <- probe_gapC_alignment_wide
probe_gapC_alignment_wide_sum$probe_gapC <- ifelse(grepl('N',probe_gapC_alignment_wide_sum$probe_gapC), 'ambiguous sequence containing "N"', probe_gapC_alignment_wide_sum$probe_gapC)

gapC_probe_plot <- ggplot(probe_gapC_alignment_wide_sum, aes(x = reorder(probe_gapC,probe_gapC,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('gapC probe')
gapC_probe_plot

# combine plots
gap_C_variant_plot <- cowplot::plot_grid(ncol = 3, gapC_F_plot, gapC_probe_plot, gapC_R_plot)


########################
### repeat same to papC#
########################

# import papC MSA
papC_alignment <- readDNAStringSet("papC.core.full.aln", format="fasta")
papC_alignment <- DNAMultipleAlignment(papC_alignment)

# locate primers
library(stringr)
papC_ref <- as.character(papC_alignment@unmasked$Reference)
#forward primer
str_locate_all(pattern ='TTTCATGGGTTGCCGGGAGTG', papC_ref)
primer_F_papC <- tidy_msa(papC_alignment, 583, 603)
#probe
str_locate_all(pattern ='TGATGCCACCTGGCTGCCTCCCT', papC_ref) #
probe_papC <- tidy_msa(papC_alignment, 672, 694)
# reverse primer
str_locate_all(pattern ='GGCATTCCCGGACTGATGCT', papC_ref) # reverse complement of primer!
primer_R_papC <- tidy_msa(papC_alignment, 709, 728)

# summarise which variant of loops occur how frequently
primer_F_papC_alignment_wide <- pivot_wider(primer_F_papC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_F_papC',2:(max(primer_F_papC$position) - min(primer_F_papC$position)+2), sep = '') 
probe_papC_alignment_wide <- pivot_wider(probe_papC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('probe_papC',2:(max(probe_papC$position) - min(probe_papC$position)+2), sep = '') 
primer_R_papC_alignment_wide <- pivot_wider(primer_R_papC, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_R_papC',2:(max(primer_R_papC$position) - min(primer_R_papC$position)+2), sep = '') 

# merge both loops
papC_alignments <- merge(primer_F_papC_alignment_wide, probe_papC_alignment_wide, by = 'name')
papC_alignments <- merge(papC_alignments, primer_R_papC_alignment_wide, by = 'name')

# check which occur how often
#remove reference
papC_alignments <- papC_alignments[papC_alignments$name != 'Reference',]
#check forward primer
count_primer_F_papC <- primer_F_papC_alignment_wide %>% 
  group_by(primer_F_papC) %>%
  summarise(n = n())
count_primer_F_papC[!grepl('N',count_primer_F_papC$primer_F_papC),] 
'TTTCATGGGTTGCCGGGAGTG'=='TTTCATGGGTTGCCGGGAGTG' #this corresponds to the primer used!
# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_F_papC_alignment_wide_sum <- primer_F_papC_alignment_wide
primer_F_papC_alignment_wide_sum$primer_F_papC <- ifelse(grepl('N',primer_F_papC_alignment_wide_sum$primer_F_papC), 'ambiguous sequence containing "N"', primer_F_papC_alignment_wide_sum$primer_F_papC)

papC_F_plot <- ggplot(primer_F_papC_alignment_wide_sum, aes(x = reorder(primer_F_papC,primer_F_papC,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('papC forward primer')
papC_F_plot


#check reverse primer
count_primer_R_papC <- primer_R_papC_alignment_wide %>% 
  group_by(primer_R_papC) %>%
  summarise(n = n())
count_primer_R_papC[!grepl('N',count_primer_R_papC$primer_R_papC),] 
'GGCATTCCCGGACTGATGCT'=='GGCATTCCCGGACTGATGCT' #this corresponds to the primer used!
# three (each 1 time, so n= 3) alternative sequences detected: GGCATTCCCGGACTGATG(T)T, GGCATTCCCGG(G)CTGATGCT, GGCATTCC(T)GGACTGATGCT  
# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_R_papC_alignment_wide_sum <- primer_R_papC_alignment_wide
primer_R_papC_alignment_wide_sum$primer_R_papC <- ifelse(grepl('N',primer_R_papC_alignment_wide_sum$primer_R_papC), 'ambiguous sequence containing "N"', primer_R_papC_alignment_wide_sum$primer_R_papC)

papC_R_plot <- ggplot(primer_R_papC_alignment_wide_sum, aes(x = reorder(primer_R_papC,primer_R_papC,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('papC reverse primer')
papC_R_plot

#check probe
count_probe_papC <- probe_papC_alignment_wide %>% 
  group_by(probe_papC) %>%
  summarise(n = n())
count_probe_papC[!grepl('N',count_probe_papC$probe_papC),] 
'TGATGCCACCTGGCTGCCTCCCT'=='TGATGCCACCTGGCTGCCTCCCT' #this corresponds to the probe used!
# 42/442 samples encode an alternative sequnece: TGATGCCACCTGGCTGCCTCC(T)T encode unambiguously a variant

# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
probe_papC_alignment_wide_sum <- probe_papC_alignment_wide
probe_papC_alignment_wide_sum$probe_papC <- ifelse(grepl('N',probe_papC_alignment_wide_sum$probe_papC), 'ambiguous sequence containing "N"', probe_papC_alignment_wide_sum$probe_papC)

papC_probe_plot <- ggplot(probe_papC_alignment_wide_sum, aes(x = reorder(probe_papC,probe_papC,function(x)-length(x)))) +
     geom_bar(stat="count", position = "dodge") +
     geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
     xlab('variant detected') +
     ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
     ggtitle('papC probe')
papC_probe_plot
# combine plots
papC_variant_plot <- cowplot::plot_grid(ncol = 3, papC_F_plot, papC_probe_plot, papC_R_plot)

##########################
### repeat same to papGII#
##########################

# import papGII MSA
papGII_alignment <- readDNAStringSet("papGII.core.full.aln", format="fasta")
papGII_alignment <- DNAMultipleAlignment(papGII_alignment)

# locate primers
library(stringr)
papGII_ref <- as.character(papGII_alignment@unmasked$Reference)
#forward primer
str_locate_all(pattern ='TCATTTCGCGAGTTACTTGG', papGII_ref)
primer_F_papGII <- tidy_msa(papGII_alignment, 543, 562)
#probe
str_locate_all(pattern ='AAGAATATCGGCGGATGCCGTC', papGII_ref) #
probe_papGII <- tidy_msa(papGII_alignment, 634, 655)
# reverse primer
str_locate_all(pattern ='CGCTAATAATCATTATGCGGC', papGII_ref) # reverse complement of primer!
primer_R_papGII <- tidy_msa(papGII_alignment, 705, 725)

# summarise which variant of loops occur how frequently
primer_F_papGII_alignment_wide <- pivot_wider(primer_F_papGII, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_F_papGII',2:(max(primer_F_papGII$position) - min(primer_F_papGII$position)+2), sep = '') 
probe_papGII_alignment_wide <- pivot_wider(probe_papGII, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('probe_papGII',2:(max(probe_papGII$position) - min(probe_papGII$position)+2), sep = '') 
primer_R_papGII_alignment_wide <- pivot_wider(primer_R_papGII, id_cols = 'name', names_from = 'position', values_from = 'character')  %>% unite('primer_R_papGII',2:(max(primer_R_papGII$position) - min(primer_R_papGII$position)+2), sep = '') 

# merge both loops
papGII_alignments <- merge(primer_F_papGII_alignment_wide, probe_papGII_alignment_wide, by = 'name')
papGII_alignments <- merge(papGII_alignments, primer_R_papGII_alignment_wide, by = 'name')

# check which occur how often
#remove reference
papGII_alignments <- papGII_alignments[papGII_alignments$name != 'Reference',]
#check forward primer
count_primer_F_papGII <- primer_F_papGII_alignment_wide %>% 
  group_by(primer_F_papGII) %>%
  summarise(n = n())
count_primer_F_papGII[!grepl('N',count_primer_F_papGII$primer_F_papGII),] 
'TCATTTCGCGAGTTACTTGG'=='TCATTTCGCGAGTTACTTGG' #this corresponds to the primer used!
# three alternative sequences were deretcted: TCATTTCGCGA(A)TTACTTGG (n=7); TCATTTCGCGA(A)TTACTT(T)G (n=3) and TCATTTCGCGAGTTACTTG(A) (n=1)
# these are  11/269 = 4.1%
# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_F_papGII_alignment_wide_sum <- primer_F_papGII_alignment_wide
primer_F_papGII_alignment_wide_sum$primer_F_papGII <- ifelse(grepl('N',primer_F_papGII_alignment_wide_sum$primer_F_papGII), 'ambiguous sequence containing "N"', primer_F_papGII_alignment_wide_sum$primer_F_papGII)

papGII_F_plot <- ggplot(primer_F_papGII_alignment_wide_sum, aes(x = reorder(primer_F_papGII,primer_F_papGII,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('papGII forward primer')
papGII_F_plot

#check reverse primer
count_primer_R_papGII <- primer_R_papGII_alignment_wide %>% 
  group_by(primer_R_papGII) %>%
  summarise(n = n())
count_primer_R_papGII[!grepl('N',count_primer_R_papGII$primer_R_papGII),] 
'CGCTAATAATCATTATGCGGC'=='CGCTAATAATCATTATGCGGC' #this corresponds to the primer used!
# detected in 392, no alternative seq detected
# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
primer_R_papGII_alignment_wide_sum <- primer_R_papGII_alignment_wide
primer_R_papGII_alignment_wide_sum$primer_R_papGII <- ifelse(grepl('N',primer_R_papGII_alignment_wide_sum$primer_R_papGII), 'ambiguous sequence containing "N"', primer_R_papGII_alignment_wide_sum$primer_R_papGII)

papGII_R_plot <- ggplot(primer_R_papGII_alignment_wide_sum, aes(x = reorder(primer_R_papGII,primer_R_papGII,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('papGII reverse primer')
papGII_R_plot


#check probe
count_probe_papGII <- probe_papGII_alignment_wide %>% 
  group_by(probe_papGII) %>%
  summarise(n = n())
count_probe_papGII[!grepl('N',count_probe_papGII$probe_papGII),] 
'AAGAATATCGGCGGATGCCGTC'=='AAGAATATCGGCGGATGCCGTC' #this corresponds to the probe used!
# an alternative sequence: AAGAATA(C)CGGCGGATGCCGTC was detected in 3/272 == 1.1% of strains
# replace all sequences with 'N' with "ambiguous sequence containing 'N'"
probe_papGII_alignment_wide_sum <- probe_papGII_alignment_wide
probe_papGII_alignment_wide_sum$probe_papGII <- ifelse(grepl('N',probe_papGII_alignment_wide_sum$probe_papGII), 'ambiguous sequence containing "N"', probe_papGII_alignment_wide_sum$probe_papGII)

papGII_probe_plot <- ggplot(probe_papGII_alignment_wide_sum, aes(x = reorder(probe_papGII,probe_papGII,function(x)-length(x)))) +
  geom_bar(stat="count", position = "dodge") +
  geom_text(stat='count', aes(label=..count.., y=..count..+10)) +
  xlab('variant detected') +
  ylab('number of isolates') + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  ggtitle('papGII probe')
papGII_probe_plot

# combine plots
papGII_variant_plot <- cowplot::plot_grid(ncol = 3, papGII_F_plot, papGII_probe_plot, papGII_R_plot)

primer_variants <- cowplot::plot_grid(gap_C_variant_plot, papC_variant_plot, papGII_variant_plot, ncol = 1)

pdf('../../../04_Presentations/PCR_variant_barplots.pdf', height = 15, width = 12)
primer_variants
dev.off()
