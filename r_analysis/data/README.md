data/
----

This directory contains raw data related to the *E.coli K-12* simluation. Usage information of all python scripts can be obtained by running the script with the `-h` option:

e.g.

```
python count_duplicated.py -h
```

# Contents:

* [count_deduplicated.py](count_duplicated.py): Given simulated reads from a genome and a ggbio Genomic Ranges file of the original CDS positions, creates a ggbio Genomic Ranges file of reads that fell within these ranges
* [ecoli/](ecoli/): *E. coli K-12 (MG1655)* (NC_000913.3) genome and coding sequence information obtained from the NCBI (March 28, 2015)
* [genomic_ranges/](genomic_ranges/): Directory containing resulting ggbio Genomic Ranges files
* [ncbi_feature_table_to_ggbio_gr.py](ncbi_feature_table_to_ggbio_gr.py): Converts an NCBI feature table file of CDS positions into a ggbio Genomic Ranges file compatible with [count_deduplicated.py](count_duplicated.py)
* [README.md](README.md): This README file
* [run.sh](run.sh): Script to run simulation experiment. Creates Genomic Ranges file from NCBI feature table ([ncbi_feature_table_to_ggbio_gr.py](ncbi_feature_table_to_ggbio_gr.py)), generates simulated *E.coli* paired-end reads [SOFAMetaSimCOG.py](SOFAMetaSimCOG.py) and creates the duplicated genomic reads file ( [count_deduplicated.py](count_duplicated.py))
* [SOFAMetaSimCOG.py](SOFAMetaSimCOG.py): Generates simluated paired-end reads from a given `.fasta` file

