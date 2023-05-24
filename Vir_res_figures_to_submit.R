# load packages 
library('ggtree')
library("phytools")
library('dplyr')
library('ggplot2')
library('GGally')
library("gridExtra")
library('stringr')
library("RColorBrewer")
library('tidyr')
library('ggnewscale')

setwd('/Users/alinecuenod/Documents/Documents_Alines_MacBook_Pro/Doc.Mobility/01_Data/03_sequences/sequence_data_outputs_to_submit')

# Import phylogenetic tree, cunstructed using RAxML
nwk <- ("RAxML_bestTree_825.raxmltree")
tree <- read.tree(nwk)
tree$tip.label <- gsub('712679-19-wh2','712679-19',  tree$tip.label)

# use mash based phylogrou assignments
phylo<-read.delim('mash_phylo.csv',quote="", sep=';',header=T, fill = TRUE)
colnames(phylo) <- c("label","best_mash_dist","Phylogroup","Phylogroup_rough") 
# rename. This strain has been resequenced, as it was relatively low coverage. 
phylo$label <- gsub('712679-19-wh2','712679-19',  phylo$label)
# save phylo all
phylo_all <- phylo

# import which sequences to consider, when only one strain per case was included
derep <- read.csv2('/Users/alinecuenod/Documents/Documents_Alines_MacBook_Pro/Doc.Mobility/01_Data/03_sequences/sequence_data_outputs_to_submit/825_strains_included.txt', sep='\t', header = F)
colnames(derep) <- 'strain'

# rename. This strain has been resequenced, as it was relatively low coverage. 
derep$strain <- gsub('712679-19-wh2','712679-19',  derep$strain)

# only look at dereplicated strains
phylo <- phylo[phylo$label %in% derep$strain,]
# extract  six digit samplenumber
phylo['TGNR']<-gsub('(\\d{6}(\\-.)*\\-\\d{2})(.*$)', '\\1\\2', phylo$label)

setdiff(tree$tip.label, phylo$label) # none are different
setdiff(phylo$label,tree$tip.label) # none are different

# re-root the tree before plotting
tree.gg<-ggtree(tree, layout = 'rectangular') 
# label nodes to choose one for re-rooting
#tree.gg + geom_text(aes(label=node), size=3)
tree2<-reroot(tree, 1301)

# plot basic treeplot
tree.gg<-ggtree(tree2, layout = 'rectangular') 
tree.gg<- (tree.gg) %<+% 
  phylo + 
  geom_tippoint(aes(col=Phylogroup), size=1, show.legend = T) + 
  theme(legend.position = c(.1,.75)) + geom_treescale(x=0.05, y= 600)
tree.gg

# Virulence and resistance factors
# use the NCBI resistance database
ncbi<-read.csv2('ncbi_resistance.tab', sep='\t')
# use the UPEC VFDB, published by Biggel et al. (https://doi.org/10.1038/s41467-020-19714-9)
vir_vf<-read.csv2('EcVGDB_virulence.tab', sep='\t')

# select virulence factors which are attributed to the class 'Invasion' or 'Ashesion / Invasion'
#vir_vf_adhe <- vir_vf[grepl('fimbr| pil|dhesion',vir_vf$PRODUCT),]
# better: check class of virulence factor
vir_vf['class'] <- gsub('(.*\\-.*class\\_)(.*$)', '\\2', vir_vf$PRODUCT)
vir_vf$class <- gsub('Invasion', 'Adhesion_invasion', vir_vf$class)

# merge resistance file and file with resistance genes
ncbi_vir<-merge(ncbi, vir_vf, by=intersect(colnames(ncbi), colnames(vir_vf)), all=T)
#filter with 95 % coverage and identity
ncbi_vir$X.COVERAGE<-as.numeric(as.character(ncbi_vir$X.COVERAGE))
ncbi_vir$X.IDENTITY<-as.numeric(as.character(ncbi_vir$X.IDENTITY))
# only include hits with a coverage and identity over 95% 
ncbi_vir<-ncbi_vir[ncbi_vir$X.COVERAGE>95 & ncbi_vir$X.IDENTITY>95,]
# add six digit samplenumber
ncbi_vir['TGNR']<-gsub("\\.f.*","",ncbi_vir$X.FILE)

# specifically check (i) AfaA-VIII, (ii) fimH, (iii) iuc, (iv) fimG and (v) csgB (curly)
ncbi_vir_spec_vir_check <-  ncbi_vir %>% 
  group_by(TGNR) %>%
  mutate(afaAVIII = any(grepl('AfaA\\-VIII', GENE))) %>%
  mutate(iuc = any(grepl('iuc', GENE))) %>%
  mutate(fimH = any(grepl('fimH', GENE))) %>%
  mutate(csgB = any(grepl('csgB', GENE))) %>% 
  select(TGNR, afaAVIII, fimH, iuc, csgB) %>%
  slice_head(n=1)

#write.csv2(ncbi_vir_spec_vir_check, '../11_abricate/vir_iuc_Afa_fimH_csgB.csv', row.names = F, quote = F)
# count how many virulence factors per class were identified
vir_class_check <-  ncbi_vir %>% 
  filter(DATABASE == "EcVGDB") %>%
  group_by(TGNR, class) %>%
  mutate(count = n())  %>%
  select(TGNR, class, count) %>%
  slice_head(n=1)

#write.csv2(vir_class_check, '../11_abricate/vir_class_check.csv', row.names = F, quote = F)

# count number of virulence and res factor per strain
n_vir<-ncbi_vir[ncbi_vir$DATABASE == 'EcVGDB',] %>% group_by(TGNR) %>% reframe(n_vir = n()) # for 5 strains no res gene was detected
n_ncbi<-ncbi_vir[ncbi_vir$DATABASE == 'ncbi',] %>% group_by(TGNR) %>% reframe(n_ncbi = n())

# check for derep strains, how often there is afaVIII
setdiff(tree$tip.label, ncbi_vir$TGNR) # there are no differences
# subset to only include one strain per clinical case
ncbi_vir <- ncbi_vir[ncbi_vir$TGNR %in% tree$tip.label,]
length(unique(ncbi_vir$X.FILE)) # these are 825 files as they should be

# merge the counting tables back together
n_ncbi_vir<-merge(n_vir, n_ncbi, by = 'TGNR', all = T)

# include only the ones included in tree (825)
length(unique(n_ncbi_vir[n_ncbi_vir$TGNR %in% tree$tip.label,]$TGNR)) # all
n_ncbi_vir<-n_ncbi_vir[n_ncbi_vir$TGNR %in% tree$tip.label,]
# merge with phylogroup
setdiff(n_ncbi_vir$TGNR, phylo$label)
n_ncbi_vir_phylo<-merge(n_ncbi_vir, phylo, by.x = 'TGNR', by.y = 'label')

# in one strain (131923-B-18), no resistance genes have been detected, set to 0
n_ncbi_vir_phylo$n_ncbi<- as.numeric(as.character(ifelse(is.na(n_ncbi_vir_phylo$n_ncbi), '0', n_ncbi_vir_phylo$n_ncbi)))

# count how many vir and res factors per phylogroup
n_phylo_<- n_ncbi_vir_phylo %>%
  group_by(Phylogroup) %>%
  reframe(quantile = scales::percent(c(0.25, 0.5, 0.75)),
            vir = quantile(n_vir, c(0.25, 0.5, 0.75)),
            res = quantile(n_ncbi, c(0.25, 0.5, 0.75)))

# plot number of vir and res factor per phylogroup
vir_plot<-ggplot(n_ncbi_vir_phylo, aes(Phylogroup, n_vir)) +
  geom_boxplot() + ylab('Number of virulence factors')
ncbi_plot <- ggplot(n_ncbi_vir_phylo, aes(Phylogroup, n_ncbi)) +
  geom_boxplot() + ylab('Number of resistance factors')

#subset the 'phylo' file, such that the tree labels correspond to the rownames
phylo_plot <- phylo[,c('label','Phylogroup')]
rownames(phylo_plot) <- phylo_plot$label
phylo_plot$label<-NULL

# summarise how many and which strains encoded papG and fimH
papG_fimH <- vir_vf %>% group_by(X.FILE) %>%
  dplyr::reframe(papG = sum(grepl('papG', GENE)), 
                   fimH = sum(grepl('fimH', GENE)))
papG_fimH <- as.data.frame(papG_fimH)
# set rownames to strains names, such that they correspond to the tree labels
rownames(papG_fimH) <- gsub('\\.fna', '', papG_fimH$X.FILE)
setdiff(rownames(phylo_plot), rownames(papG_fimH))
papG_fimH$X.FILE <- NULL
# replace '1' with the respective gene name, such that these are attributed different colors in the plot
papG_fimH$papG <- ifelse(papG_fimH$papG > 0 , 'papG', NA)
papG_fimH$fimH <- ifelse(papG_fimH$fimH > 0 , 'fimH', NA)

# check papG variants
papG_var <- vir_vf %>% group_by(X.FILE) %>%
  dplyr::reframe(papG_variant = ifelse(grepl('papG', GENE), GENE, 'no papG'))
# Some strains encode more than one papG variant
dupl_papG<- papG_var[papG_var$X.FILE %in% papG_var$X.FILE[duplicated(papG_var$X.FILE)],] #these have multiple papG variants
papG_var <- papG_var[!duplicated(paste0(papG_var$X.FILE, papG_var$papG_variant)),]
papG_var <- papG_var %>% group_by(X.FILE) %>%
  arrange(papG_variant) %>%
  mutate(papG_variant = paste0(papG_variant, collapse = "_")) %>%
  slice_head(n=1)

papG_var <- as.data.frame(papG_var)
# set labels to rownames for tree plot
rownames(papG_var) <- gsub('\\.fna', '', papG_var$X.FILE)
setdiff(rownames(phylo_plot), rownames(papG_var)) #
# merge to phylo file
papG_var_ <- papG_var
papG_var_['strain'] <- gsub('\\.fna', '', papG_var_$X.FILE)
phylo_plot_label <-phylo_plot
phylo_plot_label['label']<-rownames(phylo_plot_label)
papG_var_phylo <- merge(papG_var_, phylo_plot_label, by.x = 'strain', by.y = 'label')
# summarise how many strains encode the papG variant papGII
papG_var_phylo['papGII']<- ifelse(grepl('papGII$|papGII\\_', papG_var_phylo$papG_variant), 'papGII', 'no papGII') 

table(papG_var_phylo$papGII)
table(papG_var_phylo$papGII, papG_var_phylo$Phylogroup)

#write.csv2(papG_var_, '../06_nonribosomal-targets/papG_var.csv', row.names = F, quote = F)
papG_var$X.FILE <- NULL

# check cooccurence papG var and HdeA mass allel
hdea_mass <- read.csv('hdeA_masses.csv', sep='\t') # The mass for hdeA was calculated for 1 strain per case (n=825). for 7 of these strains, no hdeA gene was detected and no mass was calculated
colnames(hdea_mass)<-c('seq', 'mass')
hdea_mass['strain']<- gsub('(.+)(\\;.*)', '\\1', hdea_mass$seq)
setdiff(hdea_mass$strain,papG_var_$strain)

hdea_mass_pap <- merge(hdea_mass, papG_var_, by = 'strain', all = T)
# subset to one per case
setdiff(rownames(phylo_plot), hdea_mass_pap$strain) 
hdea_mass_pap <- hdea_mass_pap[hdea_mass_pap$strain %in% rownames(phylo_plot),]
# examine co-occurence of papG variants and the mass allele of hdeA
table(hdea_mass_pap$papG_variant, hdea_mass_pap$mass)
rownames(hdea_mass_pap) <- hdea_mass_pap$strain
# add phylogroup
hdea_mass_pap_phylo <- merge(hdea_mass_pap, phylo_plot, by='row.names')
rownames(hdea_mass_pap_phylo) <- hdea_mass_pap_phylo$Row.names

# summarise which strains encode iuc, which is needed for the systhesis of aerobactin
iuc <- vir_vf %>% group_by(X.FILE) %>%
  mutate(which_iuc = ifelse(any(GENE %in% c('shiF', 'iucA', 'iucB', 'iucC', 'iucD', 'iutA')),  paste(sort(unique(GENE[GENE %in% c('shiF', 'iucA', 'iucB', 'iucC', 'iucD', 'iutA')])), collapse = '_'), NA)) %>%
  group_by(X.FILE, which_iuc) %>%
  dplyr::reframe(iuc_operon_n = length(unique(GENE[GENE %in% c('shiF', 'iucA', 'iucB', 'iucC', 'iucD', 'iutA')])))

# the uic operon consists of 5 genes. if one gene in operon is missing, check which
iuc_5<-iuc[iuc$iuc_operon_n == '5',]
table(iuc_5$which_iuc)
# include only one strain per case
iuc_phylo <- iuc[gsub('\\..*$', '', iuc$X.FILE) %in% phylo$label, ]

# mlst
mlst <- read.csv2('mlst__mlst__Escherichia_coli#1__results.txt', sep='\t', header = T)
mlst <- mlst[,c('ST', 'Sample')]
mlst$Sample <- gsub('712679-19-wh2','712679-19', mlst$Sample)
# count how many ST131
st131 <- mlst[mlst$ST == "131",]
st131['TGNR']<-gsub('(\\d{6})(.*)', '\\1',st131$Sample)
# put aside original file before subsetting replacing rare ST with 'Other'
mlst_all <- mlst
# supset to include only one per case
mlst <- mlst[mlst$Sample %in% rownames(phylo_plot), ]
rownames(mlst) <- mlst$Sample
setdiff(rownames(phylo_plot), rownames(mlst))
mlst$Sample <- NULL
# count how many ST
length(unique(mlst$ST)) # 198 SR
# replace rare ST with 'other'
mlst$ST <- ifelse(mlst$ST %in% names(sort(table(mlst$ST),decreasing=TRUE)[1:8]), mlst$ST, 'Other')
mlst$ST<- factor(mlst$ST, levels = c( unique(mlst$ST[mlst$ST != 'Other']), 'Other'))

# plot tree with phylogroup
tree_p1<-gheatmap(tree.gg, phylo_plot, offset=0.00, width=0.1, 
                  legend_title="Phylogroup", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0) 
# add new color scale for next column
tree_p2.0 <- tree_p1 + new_scale_fill()
# add mlst
tree_p2<-gheatmap(tree_p2.0, mlst, offset=0.013, width=0.1, 
                  legend_title="ST", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(8, 'Set3'), 'lightgrey'), na.value = "white")
# add new color scale for next column
tree_p3.0 <- tree_p2 + new_scale_fill()
# add the number of virulence genes
n_vir <- as.data.frame(n_vir)
rownames(n_vir) <- n_vir$TGNR
n_vir$TGNR <- NULL
n_vir$n_vir <- as.numeric(n_vir$n_vir)
tree_p3<-gheatmap(tree_p3.0, n_vir, offset=0.026, width=0.1, 
                  legend_title="Number of Virulence Genes", low = "white", high = "blue", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0) 
# add new color scale for next column
tree_p4.0 <- tree_p3 + new_scale_fill()
# add the number of resistance genes
n_ncbi <- as.data.frame(n_ncbi)
rownames(n_ncbi) <- n_ncbi$TGNR
n_ncbi$TGNR <- NULL
n_ncbi$n_ncbi <- as.numeric(n_ncbi$n_ncbi)
tree_p4 <- gheatmap(tree_p4.0, n_ncbi, offset=0.039, width=0.1, 
                  legend_title="Number of Resistance Genes", low = "white", high = "orange", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0) +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p4

# export
pdf('../../../04_Presentations/cg_tree_n_factors_heatmap.pdf', width = 16, height = 10)
tree_p4 + coord_cartesian(clip = 'off') # adjust margins
dev.off()

#plot which strains include papG / fimH variant to phylogoup mlst tree
tree_p5 <- gheatmap(tree_p3.0, papG_fimH, offset=0.026, width=0.2, 
                                legend_title="", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0) + scale_fill_manual(values = c('darkorange','brown'), na.value="white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p5

# export
pdf('../../../04_Presentations/cg_tree_papfim_heatmap.pdf', width = 16, height = 10)
tree_p5 + coord_cartesian(clip = 'off') 
dev.off()

# capsule types and serotypes
serotypes <- read.csv2('EcOH__genes__EcOH__results.txt', sep='\t')
# extract o and H antigens
serotypes['O_antigen'] <- ifelse(serotypes$wzx != '-', gsub('(wzx\\-)(O[^\\_]+)(\\_.*$)', '\\2', serotypes$wzx), gsub('(wzm\\-)(O[^\\_]+)(\\_.*$)', '\\2', serotypes$wzm))
serotypes['H_antigen'] <- ifelse(serotypes$fliC != '-', gsub('(fliC\\-)(H[^\\_]+)(\\_.*$)', '\\2', serotypes$fliC), serotypes$fliC)

# subset to only one strain per clinical case
serotypes <- serotypes[serotypes$Sample %in% rownames(phylo_plot), ]
rownames(serotypes) <- serotypes$Sample
setdiff(rownames(phylo_plot), rownames(serotypes))
serotypes$Sample <- NULL
length(unique(serotypes$O_antigen)) # 98 different ones 
sum(!grepl('O', serotypes$O_antigen)) # 35 not coding any
# remplace the less frequent O-types with 'Other'
serotypes$O_antigen <- ifelse(serotypes$O_antigen == '-', NA, 
                              ifelse(serotypes$O_antigen %in% names(sort(table(serotypes$O_antigen),decreasing=TRUE)[1:8]), serotypes$O_antigen, 'Other'))
O_antigen <- serotypes %>% select(O_antigen)
O_antigen$O_antigen <- factor(O_antigen$O_antigen, levels = c( unique(O_antigen$O_antigen[O_antigen$O_antigen != 'Other']), 'Other'))

length(unique(serotypes$H_antigen)) # 36 different ones 
sum(!grepl('H', serotypes$H_antigen)) # 6 not encoding any
# replace rare H types with 'Other'
serotypes$H_antigen <- ifelse(serotypes$H_antigen == '-', NA, 
                              ifelse(serotypes$H_antigen %in% names(sort(table(serotypes$H_antigen),decreasing=TRUE)[1:8]), serotypes$H_antigen, 
                              'Other'))

H_antigen <- serotypes %>% select(H_antigen)
H_antigen$H_antigen <- factor(H_antigen$H_antigen, levels = c( unique(H_antigen$H_antigen[H_antigen$H_antigen != 'Other']), 'Other'))

# import capsule types
kaptive <- read.csv('kaptive_summary.csv')
kaptive['Sample'] <- gsub('\\.fna', '', kaptive$assembly)
# supset to only one per strain
kaptive <- kaptive[kaptive$Sample %in% rownames(phylo_plot), ]
kaptive_phylo <- phylo_plot
# add Sample variable
kaptive_phylo['Sample']<- rownames(kaptive_phylo)
# merge to phylp file
kaptive_phylo <- merge(kaptive, kaptive_phylo, by = 'Sample')
# only consider matches with > 80% coverage and no reference gene missing
kaptive_phylo$best.match <- ifelse(kaptive_phylo$best.match.cov > 80 & kaptive_phylo$no.ref.genes.missing == 0, kaptive_phylo$best.match, NA)
# count how many strains encode capsule
kaptive_phylo['binary']<- ifelse(!is.na(kaptive_phylo$best.match), 'Capsule', 'no Capsule')
table(kaptive_phylo$binary)

# summarise phylogroup to ExpeC associated or colonisation associated
kaptive_phylo['phylo_asso'] <- ifelse(grepl('B2|D|F', kaptive_phylo$Phylogroup), 'ExPEC', 
                                             ifelse(grepl('A|B1|C', kaptive_phylo$Phylogroup), 'Colonisation', 'Other'))
# count how mny of these classes encode capsule
table(kaptive_phylo$phylo_asso, kaptive_phylo$binary)
table(kaptive_phylo$Phylogroup, kaptive_phylo$binary)
# summarise which are the most frequent capsule types
sort(table(kaptive_phylo$best.match))

# replace rare capsule types with 'other'
kaptive_phylo$best.match <- ifelse(is.na(kaptive_phylo$best.match), NA,
                             ifelse(kaptive_phylo$best.match %in% names(sort(table(kaptive_phylo$best.match),decreasing=TRUE)[1:8]), kaptive_phylo$best.match, 
                               'Other'))
# subset to new file exclusively containing Capsule type for plotting
Capsule_type <- kaptive_phylo %>% select(best.match)
colnames(Capsule_type) <- 'K_type'
# summarise capsule  types
table(Capsule_type$K_type) 

# plot tree with all antigens (O- H and K- types)
# add o-antigen
tree_p6<-gheatmap(tree_p3.0, O_antigen, offset=0.026, width=0.1, 
                  legend_title="O antigen", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(7, 'Reds'), 'lightgrey'), na.value = "white")
# add new colour scale
tree_p7.0 <- tree_p6 + new_scale_fill()
# add H antigen
tree_p7 <- gheatmap(tree_p7.0, H_antigen, offset=0.039, width=0.1, 
                    legend_title="H antigen", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(8, 'Blues'), 'lightgrey'), na.value = "white")
# new colour scale
tree_p8.0 <- tree_p7 + new_scale_fill()
# add k-antigen
tree_p8 <- gheatmap(tree_p8.0, Capsule_type, offset=0.052, width=0.1, 
                    legend_title="K antigen", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(8, 'Greens'), 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p8 + theme(legend.position = 'bottom')
       
# export
pdf('../../../04_Presentations/cg_tree_antigen_heatmap.pdf', width = 16, height = 10)
tree_p8 + coord_cartesian(clip = 'off') 
dev.off()

# same tree but with legend for all
tree_p8.0 <- tree_p7 + new_scale_fill()
tree_p8b <- gheatmap(tree_p8.0, papG_fimH, offset=0.052, width=0.1, 
                    legend_title="Check H antigen", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(8, 'Greens'), 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p8b + theme(legend.position = 'bottom')


# add invasiveness 
## add which strains caused invasive infection
invas <- read.csv('pyseer_phenotype_825.tsv', sep='\t')
rownames(invas) <- invas$sample
invas$sample <- gsub('712679-19-wh2', '712679-19', invas$sample) 
setdiff(derep$strain, invas$sample)
invas <- invas[invas$sample %in% derep$strain,]
invas_p<-as.data.frame(invas[,'Invasiveness'])
colnames(invas_p)<-'Clinical.phenotype'
invas_p$Clinical.phenotype <- ifelse(invas_p$Clinical.phenotype=='1', 'Invasive UTI', 'No invasive UTI')
rownames(invas_p)<- rownames(invas)
# draw tree
tree_p9.0b <- tree_p8 + new_scale_fill()
tree_p9b <- gheatmap(tree_p9.0b , invas_p, offset=0.065, width=0.1, 
                     legend_title="Clinical phenotype", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c('darkred', 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p9b + theme(legend.position = 'bottom')
# export
pdf('../../../04_Presentations/antigen_invas.pdf', width = 15, height = 18)
tree_p9b + theme(legend.position = 'bottom')
dev.off()

# export legend
pdf('../../../04_Presentations/antigen_legend.pdf', width = 2, height = 18)
plot(cowplot::get_legend(tree_p8))
dev.off()

# draw stacked bar plot of capsule types on phylogroups
# replace NA by 'None'
kaptive_phylo$best.match <- ifelse(is.na(kaptive_phylo$best.match), 'None', kaptive_phylo$best.match)
# summarise which K- typer per phylogroup
kaptive_phylo_sum <- kaptive_phylo %>% 
  group_by(Phylogroup, best.match) %>%
  dplyr::reframe(n=n()) 
# add total count of phylogroupp
kaptive_phylo_sum_label <- kaptive_phylo_sum %>% 
  group_by(Phylogroup) %>%
  dplyr::reframe(counts=sum(n)) 
kaptive_phylo_sum_label['best.match']<-NA
# re-order
kaptive_phylo_sum$best.match <- factor(kaptive_phylo_sum$best.match, levels = c('KX03','KX29', "KX41", "KX21", "KX05", "K96", "KX31","KX02", 'Other', 'None'))

# plot
kap_phy <- ggplot(kaptive_phylo_sum, aes(fill=best.match, y=n, x=Phylogroup)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual('Capsule Types (kaptive)',values = c("#F781BF", "#4DAF4A", "#377EB8", "#984EA3", "#E41A1C","#FF7F00", "#A65628", "#FFFF33",  'white','lightgrey'), na.value = "white") +
  geom_text(data=kaptive_phylo_sum_label,aes(x = Phylogroup, y=counts, label=counts),
            position = position_dodge(0.9),vjust=-0.5)
# export
pdf('../../../04_Presentations/kaptive_per_phylo.pdf', width = 5.5, height = 3.5)
kap_phy + theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = 'top')
dev.off()

# export together with n_vir and n_res plot
pdf("../../../04_Presentations/kaptiveo_vir_res.pdf", height = 7, width = 3.5) 
grid.arrange(kap_phy + theme(legend.position = 'none'), vir_plot, ncbi_plot, ncol=1) 
dev.off() 

# add which papG variant to the tree
papG_var$papG_variant <- gsub('no papG\\_', '', papG_var$papG_variant)
pap <- papG_var

tree_p5.1 <- gheatmap(tree_p3.0, pap, offset=0.026, width=0.1, 
         legend_title="", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0) + scale_fill_manual(values = c('white', brewer.pal(6, 'RdYlGn')), na.value="white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
# new color scale
tree_p9.0 <- tree_p5.1 + new_scale_fill()
# add capsule type
tree_p9 <- gheatmap(tree_p9.0, Capsule_type, offset=0.039, width=0.1, 
                    legend_title="K antigen", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(8, 'Greens'), 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p9 
# export
pdf('../../../04_Presentations/kaptive_pap.pdf', width = 16, height = 10)
tree_p9 
dev.off()

# add which hdeA mass per papG variant
hdea_mass_pap_phylo_m<-as.data.frame(hdea_mass_pap_phylo[,'mass'])
colnames(hdea_mass_pap_phylo_m)<-'mass'
rownames(hdea_mass_pap_phylo_m)<- rownames(hdea_mass_pap_phylo)
hdea_mass_pap_phylo_m$mass <- factor(as.character(round(as.numeric(hdea_mass_pap_phylo_m$mass), 2))) 

# add hdeA mass allel
tree_p9b <- gheatmap(tree_p9.0, hdea_mass_pap_phylo_m, offset=0.039, width=0.1, 
                    legend_title="HdeA mass allele", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c(brewer.pal(3, 'Accent'), brewer.pal(3, 'BuPu'),brewer.pal(3, 'YlGn'), brewer.pal(3, 'Set2')), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p9b 

# export
pdf('../../../04_Presentations/pap_hdeA_allele.pdf', width = 16, height = 10)
tree_p9b 
dev.off()

# plot iuc
iuc_plot <- iuc[,c("X.FILE","iuc_operon_n")]
# if six iuc genes where detected, this is complete
iuc_plot$iuc_operon_n <- as.factor(ifelse(iuc_plot$iuc_operon_n == '6', 'complete', 
                                ifelse(iuc_plot$iuc_operon_n > 0 , 'incomplete', NA)))
iuc_plot <- as.data.frame(iuc_plot)
#iuc_plot$X.FILE <- gsub('103713-49', '103713-19', iuc_plot$X.FILE)
rownames(iuc_plot) <- gsub('.fna','',iuc_plot$X.FILE)
iuc_plot$X.FILE <-NULL

iuc_plot_check <- merge(iuc_plot, phylo_plot, by = 'row.names', all = F)
# count how frequently iuc was detected completely per phylogroup
table(iuc_plot_check$iuc_operon_n)
table(iuc_plot_check$Phylogroup)
table(iuc_plot_check$iuc_operon_n, iuc_plot_check$Phylogroup)

# add iuc to tree
# new color scale
tree_p10.0 <- tree_p5.1 + new_scale_fill()
# add iuc
tree_p10 <- gheatmap(tree_p10.0, iuc_plot, offset=0.039, width=0.1, 
                    legend_title="iuc operon", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c('black','darkgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p10 

pdf('../../../04_Presentations/pap_iuc.pdf', width = 5, height = 7)
tree_p10  
dev.off()

# import MIC as they were measured in routine diagnostics
MIC_wide_plot <- read.csv2('MIC_wide_plot.csv')
# interpret according to EUCAST 2021 https://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_10.0_Breakpoint_Tables.pdf
MIC_wide_plot['FOT_MIC'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Vitek.Res.GN._FOT), MIC_wide_plot$MIC_Vitek.Res.GN._FOT, NA))))
MIC_wide_plot['FOT_RES'] <- ifelse(MIC_wide_plot$FOT_MIC > 32, 'R', 'S') 
MIC_wide_plot['FOT_RES'] <- factor(MIC_wide_plot$FOT_RES, levels = c('S', 'R'))
MIC_wide_plot['NFT_MIC'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Vitek.Res.GN._NFT), MIC_wide_plot$MIC_Vitek.Res.GN._NFT, NA))))
MIC_wide_plot['NFT_RES'] <- ifelse(MIC_wide_plot$NFT_MIC > 64, 'R', 'S') 
MIC_wide_plot['NFT_RES'] <- factor(MIC_wide_plot$NFT_RES, levels = c('S', 'R'))
MIC_wide_plot['CTR_MIC_Vitek'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Vitek.Res.GN._CTR), MIC_wide_plot$MIC_Vitek.Res.GN._CTR, NA))))
MIC_wide_plot['CTR_MIC_Etest'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Etest_CTR), MIC_wide_plot$MIC_Etest_CTR, 
                                                                                  ifelse(!is.na(MIC_wide_plot$MIC_Etest.Enterobacteriaceae.URIN_CTR), MIC_wide_plot$MIC_Etest.Enterobacteriaceae.URIN_CTR, NA)))))
#use same breakpoints for etest and vitekMS MIC, all are in mg/L 
MIC_wide_plot['CTR_RES'] <- ifelse(!is.na(MIC_wide_plot$CTR_MIC_Vitek), ifelse(MIC_wide_plot$CTR_MIC_Vitek > 2,'R', ifelse(MIC_wide_plot$CTR_MIC_Vitek <= 1, 'S', 'I')), 
                                          ifelse(!is.na(MIC_wide_plot$CTR_MIC_Etest), ifelse(MIC_wide_plot$CTR_MIC_Etest> 2,'R', ifelse(MIC_wide_plot$CTR_MIC_Etest  <= 1, 'S', 'I')), NA))
                                   
MIC_wide_plot['CTR_RES'] <- factor(MIC_wide_plot$CTR_RES, levels = c('S', 'I','R'))
MIC_wide_plot['MER_MIC_Vitek'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Vitek.Res.GN._MER), MIC_wide_plot$MIC_Vitek.Res.GN._MER, NA))))
MIC_wide_plot['MER_MIC_Etest'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Etest_MER), MIC_wide_plot$MIC_Etest_MER, NA))))

MIC_wide_plot['MER_RES'] <- ifelse(!is.na(MIC_wide_plot$MER_MIC_Vitek), ifelse(MIC_wide_plot$MER_MIC_Vitek > 8,'R', ifelse(MIC_wide_plot$MER_MIC_Vitek <= 2, 'S', 'I')), 
                                          ifelse(!is.na(MIC_wide_plot$MER_MIC_Etest), ifelse(MIC_wide_plot$MER_MIC_Etest > 8,'R', ifelse(MIC_wide_plot$MER_MIC_Etest <= 2, 'S', 'I')), NA))

MIC_wide_plot['MER_RES'] <- factor(MIC_wide_plot$MER_RES, levels = c('S', 'I','R'))
MIC_wide_plot['CIP_MIC_Vitek'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Vitek.Res.GN._CIP), MIC_wide_plot$MIC_Vitek.Res.GN._CIP, NA))))
MIC_wide_plot['CIP_MIC_Etest'] <- as.numeric(as.character(gsub('=|>|<','', ifelse(!is.na(MIC_wide_plot$MIC_Etest_CIP), MIC_wide_plot$MIC_Etest_CIP, 
                                                                                  ifelse(!is.na(MIC_wide_plot$MIC_Etest.Enterobacteriaceae.URIN_CIP), MIC_wide_plot$MIC_Etest.Enterobacteriaceae.URIN_CIP, NA)))))

MIC_wide_plot['CIP_RES'] <- ifelse(!is.na(MIC_wide_plot$CIP_MIC_Vitek), ifelse(MIC_wide_plot$CIP_MIC_Vitek > 0.5,'R', ifelse(MIC_wide_plot$CIP_MIC_Vitek <= 0.25, 'S', 'I')),
                                   ifelse(!is.na(MIC_wide_plot$CIP_MIC_Etest), ifelse(MIC_wide_plot$CIP_MIC_Etest > 0.5,'R', ifelse(MIC_wide_plot$MER_MIC_Etest <= 0.25, 'S', 'I')), NA))

MIC_wide_plot['CIP_RES'] <- factor(MIC_wide_plot$CIP_RES, levels = c('S', 'I','R'))

MIC_wide_plot['label'] <- paste0(MIC_wide_plot$TGNR, '-', gsub('^20', '', MIC_wide_plot$year))
MIC_wide_plot <- as.data.frame(MIC_wide_plot[MIC_wide_plot$label %in% rownames(pap),c("label", "CTR_RES", "MER_RES", "FOT_RES", "NFT_RES","CIP_RES")])
MIC_wide_plot <- MIC_wide_plot[!duplicated(MIC_wide_plot),]

#for these multiple AMR profiles were acquired, remove, as these cannot be unambiguously be attributed to a strain sequenced for this study
MIC_wide_plot <- MIC_wide_plot[!MIC_wide_plot$label %in% c('701336-20', '715392-19', '719670-19', '719819-18', '720902-18', '721003-18'),]
rownames(MIC_wide_plot) <- MIC_wide_plot$label
MIC_wide_plot$label <- NULL

MIC_wide_mlst <- merge(MIC_wide_plot, mlst, by = "row.names", all.y = T)

tree_p11 <- gheatmap(tree_p3.0, MIC_wide_plot, offset=0.026, width=0.40, 
                    legend_title="Phenotypic AMR", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c('Yellow', 'Red', 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p11 

# export
pdf('../../../04_Presentations/amr_phenotypic.pdf', width = 9, height = 7)
tree_p11  
dev.off()

# subset only ceftriaxone resistance
MIC_wide_plot_ctr<-as.data.frame(MIC_wide_plot[,'CTR_RES'])
colnames(MIC_wide_plot_ctr) <- "Ceftriaxone resistance"
rownames(MIC_wide_plot_ctr)<-rownames(MIC_wide_plot)
MIC_wide_plot_ctr_mlst <- merge(MIC_wide_plot_ctr, mlst, by = "row.names", all.y = T)
# add to mlst
MIC_wide_plot_mlst <- merge(MIC_wide_plot, mlst, by = "row.names", all.y = T)


#plot papG variant and iuc in one plot with tree
tree_p10b <- tree_p10 + new_scale_fill()
tree_p10b <- gheatmap(tree_p10b, MIC_wide_plot_ctr, offset=0.052, width=0.1, 
                     legend_title="Ceftriaxone resistance", colnames = TRUE, colnames_position = 'top', colnames_angle = 60, font.size = 6, hjust = 0)  + scale_fill_manual(values = c('Red', 'lightgrey'), na.value = "white") +
  geom_treescale(width = 0.1, label = 'Substitutions per site', offset.label = 15)
tree_p10b 

pdf('../../../04_Presentations/pap_iuc_ctr.pdf', width = 6, height = 7)
tree_p10b   
dev.off()

all(invas$sample %in% papG_var_phylo$strain)
invas_pap <- merge(invas, papG_var_phylo, by.x = 'sample', by.y = 'strain')
table(invas_pap$papGII, invas_pap$Invasiveness)
# calculate specificity, sensitivity, PPR and NPR when using papGII as single predictor for invasiveness
# sensitivity: TP / (TP + FN)
97/(97+164)
# specificity: TN / (TN + FP)
490/(490+74)
# positive predictive value: TP / (TP+FP)
97/(97+74)
# negative predictive value: TN / (TN+FN)
490/(490+164)
# Accuracy: TP + TN / (TP+ FP + TN + FN) 
(97+490)/825


# repeat the same exclusively for matching samples (same strain in urine than in bloodstream)
derep_strict <- read.csv2('pyseer_phenotyper_only_matching.tsv', sep = '\t')
setdiff(derep_strict$sample, invas_pap$sample)
invas_pap_strict <- invas_pap[invas_pap$sample %in% derep_strict$sample,]
table(invas_pap_strict$papGII, invas_pap_strict$Invasiveness)
# calculate specificity, sensitivity, PPR and NPR when using papGII as single predictor for invasiveness
# sensitivity: TP / (TP + FN)
52/(52+41)
# specificity: TN / (TN + FP)
490/(490+74)
# positive predictive value: TP / (TP+FP)
52/(52+74)
# negative predictive value: TN / (TN+FN)
490/(490+41)
# Accuracy: TP + TN / (TP+ FP + TN + FN) 
(52+490)/657


