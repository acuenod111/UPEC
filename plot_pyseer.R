library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(ggplot2)

# adapted from https://pyseer.readthedocs.io/en/master/tutorial.html#interpreting-significant-k-mers

# import gene hits, annoated by pyseer
gene_hits_sign_825 = read.table("gene_hits_sign_825.txt", stringsAsFactors=FALSE, header=TRUE)
gene_hits_sign_657 = read.table("gene_hits_sign_657_update.txt", stringsAsFactors=FALSE, header=TRUE)


# plot significant hits
gene_hits_sign_825['gene_label'] <-  gene_hits_sign_825$gene
p1_825<-ggplot(gene_hits_sign_825, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive Infection, all cases (n=825)") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)") + theme(legend.position = "bottom")
#p1_825

gene_hits_sign_657['gene_label'] <- ifelse(gene_hits_sign_657$maxp> 12, gene_hits_sign_657$gene, '')
p1_657<-ggplot(gene_hits_sign_657, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene_label)) +
  geom_point(alpha=0.5) +
  geom_text_repel(aes(size=60), show.legend = FALSE, colour='black') +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive Infection, UT as port of entry (n=657)") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)") + theme(legend.position = "bottom")
#p1_657


gene_hits_sign_825['GWAS']<-'All Cases (n=825)'
gene_hits_sign_657['GWAS']<-'Urinary Tract as Port of Entry (n=657)'

gene_hits_sign_both <- merge(gene_hits_sign_657, gene_hits_sign_825, by = colnames(gene_hits_sign_657), all = T)
gene_hits_sign_both$gene_label <- gsub('100033-19', 'A', gene_hits_sign_both$gene_label)
gene_hits_sign_both$gene_label <- gsub('100888-20', 'B', gene_hits_sign_both$gene_label)

p1_both<-ggplot(gene_hits_sign_both, aes(x=avg_beta, y=maxp, colour=avg_maf, size=hits, label=gene_label)) +
  geom_point(alpha=0.5) +
  facet_wrap(. ~ GWAS, scales='free') +
  geom_text_repel(aes(size=30), show.legend = FALSE, colour='black') +
  scale_size("Number of k-mers", range=c(1,10)) +
  scale_colour_gradient('Average MAF') +
  theme_bw(base_size=14) +
  ggtitle("Invasive Infection") +
  xlab("Average effect size") +
  ylab("Maximum -log10(p-value)") + theme(legend.position = "bottom")
p1_both

pdf('../../../04_Presentations/pyseer_sign_825_657_update.pdf', height=6, width = 8)
p1_both
dev.off()
