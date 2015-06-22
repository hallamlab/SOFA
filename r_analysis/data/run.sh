#!/bin/bash

# create genomic ranges from E.coli feature table
python ncbi_feature_table_to_ggbio_gr.py -i ecoli/NC_000913.3.txt -o genomic_ranges/NC_000913.3_CDS_gr.txt

# create simulated E.coli paired-end reads
python SOFAMetaSimCOG.py -i ecoli/NC_000913.3.fasta -o ecoli/EcoliSim2.fasta -n 10000 -l 150 -t 301 -r 0

# check for duplicates
python count_duplicated.py -r ecoli/EcoliSim2.fasta -g genomic_ranges/NC_000913.3_CDS_gr.txt -m 30 -o genomic_ranges/EcoliSim2_dups_gr.txt
