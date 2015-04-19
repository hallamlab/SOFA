# SOFA: Short-ORF Functional Annotation Pipeline

Aria S. Hahn, Niels W. Hanson, Dongjae Kim, Kishori M. Konwar and Steven J. Hallam

Accurate description of the microbial communities driving matter and energy transformation in complex ecosystems such as soils cannot yet be effectively accomplished using assembly-based approaches despite the rise of next generation sequencing technologies. Here we present SOFA, a modular open source pipeline enabling comparative functional annotation of unassembled short-read data.  The pipeline attempts to merge mate pairs in fastq files, predicts open reading frames (ORFs) on merged and unmerged reads as small as 70 bps, and completes an additional step, we term 'deduplication'. Deduplication prevents the double counting of ORFs predicted from unmerged paired-end reads by checking for homologous annotations that span the same gene, allowing for quantitatively accurate gene counts. SOFA enables downstream processing stages within the existing [MetaPathways pipeline](https://github.com/hallamlab/MetaPathways2), and is available for download as a stand alone application at [https://github.com/hallamlab/SOFA](https://github.com/hallamlab/SOFA).

More information can be found on the [wiki](https://github.com/hallamlab/SOFA/wiki)!

Source code for the binaries required by SOFA are included. Please run make on the project root directly to generate the required binaries.

Binary: Website: Repository      
FLASH: http://ccb.jhu.edu/software/FLASH/ : http://sourceforge.net/projects/flashpage/files/
LAST: http://last.cbrc.jp/ : https://github.com/dongjaek/LAST
FragGeneScanPlus: http://hallam.microbiology.ubc.ca/ : https://github.com/dongjaek/FragGeneScan

