require(ggplot2)
require(ggrepel)
library(ggrepel)
library(cowplot)
library(dplyr)
library(ggplot2)

# adapted from https://pyseer.readthedocs.io/en/master/tutorial.html#interpreting-significant-k-mers

# import gene hits, annoated by pyseer
gene_hits = read.table("/Users/aline/Doc.Mobility/01_Data/03_sequences/07_pyseer/gene_hits_all_828_no_ref.txt", stringsAsFactors=FALSE, header=TRUE)

# 
#gene_hits_5 <- gene_hits[gene_hits$maxp >= 4,]
#write.table(gene_hits_5, '/Users/aline/Doc.Mobility/01_Data/03_sequences/07_pyseer/gene_hits_4_no_ref.csv', row.names = F, dec = '.', quote = F, sep = ';')

# plot significant hits
p1<-ggplot(gene_hits_5, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive Infection") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)") + theme(legend.position = "bottom")

# label only the ones with a logpvalue > 4.3
gene_hits['gene_label'] <- ifelse(gene_hits$maxp> 4.3, gene_hits$gene, '')

# plot all hits
p1b<-ggplot(gene_hits, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene_label)) +
  geom_point(alpha=0.5) +
  scale_size("Number of k-mers", range=c(1,10)) +
  geom_text_repel(aes(size=60), show.legend = FALSE, colour='black', max.overlaps = 30) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive Infection") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)") + theme(legend.position = "bottom") +
  geom_hline(yintercept = -log10(2.78E-08), col = "red", linetype = "dashed")

pdf('/Users/aline/Doc.Mobility/04_Presentations/pyseer_828.pdf', height=6, width = 8)
p1b
dev.off()
