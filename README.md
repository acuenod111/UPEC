This directory contains scripts which were used to analyse bacterial data included in the study "Bacterial genome wide association study substantiates papGII of E. coli as a patient independent driver of urosepsis"

# Vir_res_figures.R
This script combines the output of different tools which were used to characterise E. coli strains (all files can be found in 'sequence_data_output' ). This includes: 
- mash_phylo.csv: Assignment of each strain to one out of 14 phylogroups. Phylogroups were assigned by calculating Mash distances to 14 reference strains (adapted from https://doi.org/10.1038/s42003-020-01626-5)
- mlst__mlst__Escherichia_coli#1__results: Multi-Locus Sequence Types as assigned by srst2 (https://github.com/katholt/srst2)
- EcOH__genes__EcOH__results.txt: O- and H- antigens as assigned by srst2 (https://github.com/katholt/srst2)
- kaptive_summary.csv: Capsule types as assigned by 'fastKaptive' (https://github.com/rmostowy/fastKaptive) 
- hdeA_masses.csv: The predicted mass of the protein HdeA using transseq (https://www.expasy.org/resources/protparam)
- ncbi_resistance.tab: resistance genes identified  using the NCBI AMRFinder database (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047) via abricate (https://github.com/tseemann/abricate) 
- EcVGDB_virulence.tab: virulence genes identified using the EcVGDB database (https://github.com/MBiggel/UPEC_study) via abricate (https://github.com/tseemann/abricate) 
- MIC_wide_plot.csv: Minimal inhibitory concentrations for Ceftriaxone, Fosfomycine. Nitrofurantoine, Meropenem and Ciprofloxacine, measured in routine diagnostics using micro-dilution assays. 
- RAxML_bestTree_825.raxmltree: A phylogenetic tree constructed using RAxML and including one strain per clinical case (n=825)

The file '825_strains_included.txt' lists one representative strain per clinical case. 

# plot_pyseer.R
This script was adapted from https://pyseer.readthedocs.io/en/master/tutorial.html#interpreting-significant-k-mers and plots the association and average effect size of all unitigs with the endpoint 'invassive infection'. 
One input file is required: 
- gene_hits_all_825_no_ref.txt: gene hits identified by pyseer (https://pyseer.readthedocs.io/en/master/)

# specta_preprocessing-bruker.R / specta_preprocessing-shimadzu.R
These files pick peaks from the raw spectra, which were either generated using a mass spectrometer from Shimadzu (mzXML files) or from Bruker (fid files)
These scripts can be run with

Rscript specta_preprocessing-shimadzu.R ./mzXml-processed_Launchpad ./Shimadzu/csv_median ./poso.tgnr.invas.csv
Rscript specta_preprocessing-bruker.R ./fid ./Bruker/csv_median

where /fid and /mzXml-processed_Launchpad indicate the directory to the rawfiles, /csv_median the path to the output files and "poso.tgnr.invas.csv" a file translating between target plate positions and samplenames required for the Shimadzu pre-processing. 
The peak picking is based on the packages MALDIQuantForeign and MALDIQuant (https://github.com/sgibb/MALDIquant)

# binary_from_spectra.R
This script plots the presence / absence of MALDI-TOF MS peaks summarised by phylogroup and papG variant. 
It requires the following input files: 
- one_per_case_strain_825.txt: A list of all paths to the 828 assemblies
- pagG_var.csv: papG variants as assigned by the EcVGDB virulence database. 
- mash_phylo.csv: Assignment of each strain to one out of 14 phylogroups. Phylogroups were assigned by calculating Mash distances to 14 reference strains (adapted from https://doi.org/10.1038/s42003-020-01626-5)

# PCR_snippy_eval.R
This script evaluates which variants where detected for PCR primers and probes, using the variant caller Freebayes via snippy (https://github.com/tseemann/snippy)
It requires the full alignments which where outputted by snippy as input files:
- gapC.core.full.aln
- papC.core.full.aln
- papGII.core.full.aln

