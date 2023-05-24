#!/usr/bin/env python
# coding: utf-8

import Bio
from Bio.SeqUtils import molecular_weight
from Bio import SeqIO
import csv
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

input_file = "hdeA_aa_nostop.fas"
hdeA = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

# hdeA has a 21 aa long signal peptide, cut this
masses = dict()
# calculate mass
for header, value in hdeA.items():
    masses[header] =  molecular_weight(hdeA[header].seq[21:], "protein")
# export
massesdf = pd.DataFrame.from_dict(masses,orient='index')
massesdf.to_csv( "hdeA_masses.csv", sep='\t', encoding='utf-8')


