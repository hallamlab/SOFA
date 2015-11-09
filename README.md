# SOFA: Short-ORF Functional Annotation Pipeline

Aria S. Hahn, Niels W. Hanson, Dongjae Kim, Kishori M. Konwar and Steven J. Hallam

Accurate description of the microbial communities driving matter and energy transformation in complex ecosystems such as soils cannot yet be effectively accomplished using assembly-based approaches despite the rise of next generation sequencing technologies. Here we present SOFA, a modular open source pipeline enabling comparative functional annotation of unassembled short-read data.  The pipeline attempts to merge mate pairs in fastq files, predicts open reading frames (ORFs) on merged and unmerged reads as small as 70 bps, and completes an additional step, we term 'deduplication'. Deduplication prevents the double counting of ORFs predicted from unmerged paired-end reads by checking for homologous annotations that span the same gene, allowing for quantitatively accurate gene counts. SOFA enables downstream processing stages within the existing [MetaPathways pipeline](https://github.com/hallamlab/MetaPathways2), and is available for download as a stand alone application at [https://github.com/hallamlab/SOFA](https://github.com/hallamlab/SOFA).

More usage information can be found on the [wiki](https://github.com/hallamlab/SOFA/wiki)!

See [r_analysis/](r_analysis/) for more information on downstream analysis.

If using SOFA for academic work, please cite:

Aria S. Hahn, Niels W. Hanson, Dongjae Kim, Kishori M. Konwar, Steven J. Hallam. *Assembly independent functional annotation of short-read data using SOFA: Short-ORF Functional Annotation*, Proceedings of the 2015 IEEE Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2015), Niagara Falls, Canada, August 12-15, 2015. doi:[10.1109/CIBCB.2015.7300324](http://dx.doi.org/10.1109/CIBCB.2015.7300324).
