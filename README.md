[TOC]

# EasyNanoMeta: A pipeline for nanopore long-read metagenomic analysis

With the widespread application of nanopore sequencing in microbiology, there has been a significant increase in the generation of long-read sequencing data. This poses a challenge in effectively mining valuable biological information from the vast amount of nanopore sequencing data using bioinformatics tools. Here, we systematically summarizes bioinformatics analysis methods for nanopore long-read metagenomics, validating the practical implementation of these workflows using various simulated and experimental datasets from a one health perspective. In addition, we comprehensively evaluated the performance and computational requirements of various bioinformatic tools.

## Workflow overview

![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure1.jpg)

This comprehensive statement encompasses the various steps involved in the analysis process, along with the corresponding tools employed. The initial step in analyzing nanopore sequencing data involves data quality control, including basecalling, filtering and host sequences removal. Then, the the main steps of nanopore metagenomic data analysis were summarized into assembly-free strategy and assembly-based strategy. Assembly-free long-read metagenomic profiling can quickly estimate the abundance of species and functional genes employing few computing resources. In comparison, assembly-based nanopore metagenomic data analysis enable us to study more high-quality MAGs and their functions.

## How to use?  
### 1. Use the pipeline through Singularity
· At first, following the document (https://github.com/sylabs/singularity/blob/main/INSTALL.md) to install Singularity.  
· Then, download the EasyNanoMeta.sif in your computer (https://figshare.com/articles/software/A_singularity_sandbox_for_EasyNanoMeta_/27014869?file=49175110).  
· Using script easynanometa.py to perform data analysis.
```
./easynanometa.py -h
usage: easynanometa.py [-h] -f FOLDER [-t THREADS] -host-removal-reference HOST_REMOVAL_REFERENCE
                       -centrifuge-db CENTRIFUGE_DB -kraken2-db KRAKEN2_DB
                       [-checkm2-db CHECKM2_DATABASE] [-gtdbtk-db GTDBTK_DATABASE]

Execute the EasyNanoMeta-pipeline for nanopore metagenomic data analysis.

options:
  -h, --help            show this help message and exit
  -f FOLDER, --folder FOLDER
                        Absolute path of the folder to search for fastq/fq files.
  -t THREADS, --threads THREADS
                        Number of threads to use for all operations (default: 24).
  -host-removal-reference HOST_REMOVAL_REFERENCE, --host-removal-reference HOST_REMOVAL_REFERENCE
                        Path to the reference host genome fasta file for host removal.
  -centrifuge-db CENTRIFUGE_DB, --centrifuge-db CENTRIFUGE_DB
                        Path to centrifuge database path.
  -kraken2-db KRAKEN2_DB, --kraken2-db KRAKEN2_DB
                        Path to Kraken2 database path.
  -checkm2-db CHECKM2_DATABASE, --checkm2-database CHECKM2_DATABASE
                        Path to checkm2 database path.
  -gtdbtk-db GTDBTK_DATABASE, --gtdbtk-database GTDBTK_DATABASE
                        Path to GTDB-Tk database path.
```

Configure databases for easynanometa.sif
```
mkdir ~/db
cd ~/db

#human_genome:
wget https://ftp.ncbi.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/C
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz

#kraken2:
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
tar -zcvf k2_standard_20230605.tar.gz -C ~/db/k2_standard/

#centrifuge:
wget https://zenodo.org/record/3732127/files/h%2Bp%2Bv%2Bc.tar.gz?download=1
tar -zxvf centrifuge_h+p+v.tar.gz

#checkm2
wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz
tar -zxvf checkm2_database.tar.gz

#gtdbtk
wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
tar -zxvf gtdbtk_data.tar.gz
```
Put ``easynanometa.sif`` and ``easynanometa.py`` in the same floder.   
Example usage:
```
easynanometa.py -f /path/to/folder -t 40 \
-host-removal-reference /path/to/reference_genome.fasta \
-centrifuge-db /path/to/database/centrifuge/hpvc/hpvc \
-kraken2-db /path/to/database/kraken2_db/k2-standard \
-checkm2-db /path/to/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd \
-gtdbtk-db /path/to/database/gtdbtk/release214
```

Example input folders: human_sputums
```
mkdir human_sputums
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR008/8641/SRR8641382/SRR8641382.lite.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR008/8641/SRR8641382/SRR8641383.lite.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR008/8641/SRR8641382/SRR8641384.lite.1
fastq-dump SRR864138*.lite.1
```

Example output folders: easynanometa_result
```
easynanometa_result
├── abricate_out
│   ├── abricate_arg_out
│   ├── abricate_is_out
│   ├── abricate_vf_out
│   ├── arg_abundance_out
│   ├── SRR8641382.lite.1.fasta
│   ├── SRR8641383.lite.1.fasta
│   └── SRR8641384.lite.1.fasta
├── adapters_removal_out
│   ├── SRR8641382.lite.1_output.fastq
│   ├── SRR8641383.lite.1_output.fastq
│   └── SRR8641384.lite.1_output.fastq
├── centrifuge_out
│   ├── SRR8641382.lite.1_fitted_raw_kraken_report
│   ├── SRR8641382.lite.1_report
│   ├── SRR8641382.lite.1_result
│   ├── SRR8641383.lite.1_fitted_raw_kraken_report
│   ├── SRR8641383.lite.1_report
│   ├── SRR8641383.lite.1_result
│   ├── SRR8641384.lite.1_fitted_raw_kraken_report
│   ├── SRR8641384.lite.1_report
│   └── SRR8641384.lite.1_result
├── checkm2_out
│   └── SRR8641382.lite.1_checkm2_out
├── gtdbtk_out
│   └── SRR8641382.lite.1_gtdbtk_out
├── host_removal_out
│   ├── human_genome.min
│   ├── SRR8641382.lite.1.fasta
│   ├── SRR8641382.lite.1_fitted_raw.fastq
│   ├── SRR8641382.lite.1_unique.fastq
│   ├── SRR8641383.lite.1.fasta
│   ├── SRR8641383.lite.1_fitted_raw.fastq
│   ├── SRR8641383.lite.1_unique.fastq
│   ├── SRR8641384.lite.1.fasta
│   ├── SRR8641384.lite.1_fitted_raw.fastq
│   └── SRR8641384.lite.1_unique.fastq
├── kraken2_out
│   ├── SRR8641382.lite.1_kraken2_report
│   ├── SRR8641382.lite.1_kraken2_result
│   ├── SRR8641383.lite.1_kraken2_report
│   ├── SRR8641383.lite.1_kraken2_result
│   ├── SRR8641384.lite.1_kraken2_report
│   └── SRR8641384.lite.1_kraken2_result
├── metaflye_out
│   ├── SRR8641382.lite.1_flye_out
│   ├── SRR8641383.lite.1_flye_out
│   └── SRR8641384.lite.1_flye_out
├── nextpolish_out
│   ├── SRR8641382.lite.1.cfg
│   ├── SRR8641382.lite.1.fofn
│   ├── SRR8641382.lite.1_nextpolish_out
│   ├── SRR8641383.lite.1.cfg
│   ├── SRR8641383.lite.1.fofn
│   ├── SRR8641383.lite.1_nextpolish_out
│   ├── SRR8641384.lite.1.cfg
│   ├── SRR8641384.lite.1.fofn
│   └── SRR8641384.lite.1_nextpolish_out
└── semi_bin_out
    ├── bam_out
    └── SRR8641382.lite.1_bin_out
```

Alternatively, using ``easynanometa2.py``to call individual tool.
```
./easynanometa2.py -h
usage: easynanometa2.py [-h]
                {flye,kraken2,abricate,adapters-removal,host-removal,centrifuge,arg-abundance,nextpolish,semibin,checkm2,gtdbtk}
                ...

Batch execute Flye, Kraken2 or abricate for multiple fastq files.

positional arguments:
  {flye,kraken2,abricate,adapters-removal,host-removal,centrifuge,arg-abundance,nextpolish,semibin,checkm2,gtdbtk}
                        Choose the tool to run.
    flye                Run Flye assembler
    kraken2             Run Kraken2 classifier
    abricate            Run abricate to identify functional genes.
    adapters-removal    Run porechop_abi to remove adapters.
    host-removal        Run host_removal to remove host genome.
    centrifuge          Run Centrifuge classifier.
    arg-abundance       Run calculating abundance for ARGs.
    nextpolish          Run polish for flye results.
    semibin             Run SemiBin for NextPolish results.
    checkm2             Run Checkm2 for SemiBin output.
    gtdbtk              Run GTDBTK for SemiBin output.

options:
  -h, --help            show this help message and exit
```

### 2. Use the pipeline through shell scripts
#### Readme.md  
Introduction  
#### install.sh  
Dependencies installation. Folowing the script to install tools for nanopore metagenomic analysis.  
#### pipeline.sh  
Command-line analysis for Linux. A guidance for analyzing nanopore metagenomic data.  
#### Python scripts for data analysis  
We provided many scripts for data processing.  
For examples:  
Functional gene abundance calculation:  
``$ python abundance_calculate.py --help``  
``usage: abundance_calculate.py [-h] [--i I] [--data_size DATA_SIZE] [--title TITLE] [--p P] [--output OUTPUT]``  

``options:``  
``  -h, --help            show this help message and exit``  
``  --i I, -i I           Input data.``  
``  --data_size DATA_SIZE, -d DATA_SIZE``  
``                        --data_size, -d``  
``  --title TITLE``  
``  --p P, -p P           The prefix of result.``  
``  --output OUTPUT, -o OUTPUT``  
``                        Output direction.``  

Co-located analysis:  
``$ python co-located.py --help``  
``usage: co-located.py [-h] [--i I]``  
  
``options:``  
``  -h, --help  show this help message and exit``  
``  --i I``  

Merge multiple tables:  
``$ python merge_tables.py -h``  
``usage: merge_tables.py [-h] [--input_dir INPUT_DIR] [--identifier IDENTIFIER]``  
``                       [--prefix PREFIX] [--output OUTPUT] [--column COLUMN]``  

``optional arguments:``  
``  -h, --help            show this help message and exit``  
``  --input_dir INPUT_DIR, -i INPUT_DIR``  
``                        Input data directory.``  
``  --identifier IDENTIFIER, -d IDENTIFIER``  
``                        The uniform identifier for name of the prepared files.``  
``  --prefix PREFIX, -p PREFIX``  
``                        The prefix of output result.``  
``  --output OUTPUT, -o OUTPUT``  
``                        Output data directory.``  
``  --column COLUMN, -c COLUMN``  
``                        Column name used as the key for merging between``  
``                        Dataframes.``  


#### Rscripts_for_ploting  
Examples for ploting using R.    

## Dependencies
### Software management
**Miniconda3** (https://docs.anaconda.com/free/miniconda/index.html)  
### Host removal and quality control
**Minimap2** (https://github.com/lh3/minimap2)  
**samtools** (https://github.com/samtools/samtools)  
**bedtools** (https://github.com/arq5x/bedtools2)  
**fastp** (https://github.com/OpenGene/fastp)  
### Taxonomic and functional composition identification
**Centrifuge** (http://www.ccb.jhu.edu/software/centrifuge)  
**Kraken2** (https://github.com/DerrickWood/kraken2)  
**Abricate** (https://github.com/tseemann/abricate)  
### Long-read metagenomic assembly
**MetaFlye** (https://github.com/fenderglass/Flye)  
**Canu** (https://github.com/marbl/canu)  
**wtdbg2** (https://github.com/ruanjue/wtdbg2)  
**NextDenove** (https://github.com/Nextomics/NextDenovo)  
### Hybrid metagenomic assembly
**OPERA-MS** (https://github.com/CSB5/OPERA-MS)  
**MetaSPAdes** (https://github.com/ablab/spades)  
**MetaPlatanus** (https://github.com/rkajitani/MetaPlatanus)  
**Unicycler** (https://github.com/rrwick/Unicycler)  
### Polishing  
**NextPolish** (https://github.com/Nextomics/NextPolish)  
**Pilon** (https://github.com/broadinstitute/pilon)  
### Binning
**bwa** (https://github.com/lh3/bwa)  
**Semibin** (https://github.com/BigDataBiology/SemiBin)  
**vamb** (https://github.com/RasmussenLab/vamb)  
**MetaWrap** (https://github.com/bxlab/metaWRAP)  
### Annotation
**Checkm** (https://github.com/Ecogenomics/CheckM)  
**Checkm2** (https://github.com/chklovski/CheckM2)  
**Gtdbtk** (https://ecogenomics.github.io/GTDBTk/index.html)  
**Prokka** (https://github.com/tseemann/prokka)  

## Practical guidances
### 1. A practical guidance for assembly-free long-read metagenomic data analysis
**Figure1**
![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure2.jpg)

**Statistical and visualization analysis for assembly-free long-read metagenomic data**  
We used a set of datasets from TDT human gut to perform assembly-free analysis. After obtaining the taxonomic and resistance abundance tables, we can perform many statistical and visualization analyses, including composition, diversity, and correlation. (A) The workflow for assembly-free nanopore metagenomic data analysis. (B) The microbial compositions at the phylum level. An increase of Bacteroidetes was observed in Severe group. (C) The distribution of the top 20 abundant ARGs in different samples. ARGs in Severe samples are more abundant than in Mild samples. (D) Alpha diversity analysis for microbial and ARG compositions. Alpha diversity showed that Severe group has more abundant ARGs. (E) Beta diversity analysis for ARG compositions. PCoA based on Bray-Curtis showed that there is a significant difference in ARG compositions between Mild and Severe. (F) The correlation analysis of microbial and ARG compositions. In the heat map, dark blue color blocks indicate positive correlations between different ARGs based on the Pearson Correlation Coefficient. The attachment of ARGs and the bacterial phylum indicate the correlation. (G) The correlation analysis of different ARGs based on their genetic context. The co-occurrence network diagram shows the correlation of different ARGs based on their genetic location. The nodes represent different ARGs. Connections between nodes indicate that they are related. The edge widths represent the correlation coefficient.

### 2. A practical guidance for assembly-based long-read metagenomic data analysis
**Figure2**
![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure3.jpg)

**Evaluation of metagenome assembly and binning results and computing resource consumption corresponding to each tested tool for different datasets**   
Long-read metagenomic assembly was performed with MetaFlye, Canu, NextDenove and wtdbg2. Hybrid metagenomic assembly was performed with MetaSPAdes, OPERS-MS and MetaPlatanus. Two short-read metagenomic assemblers (MetaSPAdes and MEGAHIT) are optional for OPERS-MS to carry out short-read pre-assembly. Five different datasets were evaluated. Contigs shorter than 1 kbp are not included in the statistics. Four assembly tools and two binning tools were chosen to perform assembly and binning. According to completeness and contamination, MAGs are classified as high-quality (more than 90% completeness and less than 5% contamination) and medium-quality (more than 50% completeness and less than 10% contamination). (A) The workflow for assembly-based nanopore metagenomic data analysis. (B,C) The contigs count and assembly size generated by different assembly strategies for the five datasets. For long-read assembly strategies, MetaFlye always produces the biggest assembly size, followed by wtdbg2. The assembly size does not correspond to the number of contigs, indicating that the contig length distribution is various in different assemblies. For hybrid assembly strategies, hybrid assembly based on SPAdes usually produces a bigger assembly size and a higher contigs count (MetaSPAdes and OPERA-MS_SPAdes). (D) Proportion of the assembly size with different lengths of contigs. Long-read assembly produces longer contigs than hybrid assembly. Although NextDenove usually produces contigs with a length larger than 50 kbp, the assembly size was significantly smaller than that of other assemblers. Hence, NextDenove is not suitable for metagenomic assembly. For hybrid assembly strategies, MetaPlatanus produces the most contiguous assemblies. Compared with MetaSPAdes, OPERS-MS could significantly improve the length of contigs. (F,F) The computational requirements of different assembly strategies. We used 24 cpus to test different assembly tools. Compared with the other three long-read assemblers, Canu consumed several times of computer time. MetaFlye and wtdbg2 consumed less memory than Canu and NextDenove. For hybrid assembly strategies, OPERA-MS based on MEGAHIT consumed the least amount of computer time and memory. (G,H,I) The influence of long-read only assembly and binning methods on the quality of MAGs. We take the human gut dataset as an example. The number of MAGs from MetaFlye assemblies is the highest. SemiBin shows better performance in long-read metagenomic binning than vamb. The mean contigs length of MAGs derived from NextDenove assemblies is significantly longer than that of other MAGs. In the words, NextDenove could produce more contiguous MAGs. The contigs N50 of MAGs derived from NextDenove assemblies are significantly larger than other MAGs. Overview, the combination of MetaFlye and SemiBin could produce more MAGs. Although NextDenove and SemiBin could produce longer contigs of MAGs, they lost too many MAGs.

**Figure3**
![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure4.jpg)

**The influence of different assembly strategies on the accuracy of MAGs**   
We used four hybrid assembly strategies and four long-read assembly strategies to perform metagenomic assembly of a Mock2 dataset. The long-read assemblies were polished with long-read data and combined with long-read and short-read data. The suffix “lp” indicates long-read polishing. The suffix “slp” indicates short-read and long-read polishing. After that, we extracted MAGs from different assemblies and compared MAGs with reference genomes. (A) Workflow for assessing the quality of MAGs. We compared the MAGs recovered using different analysis strategies with reference genomes to assess the quality of MAGs. (B) The genome size difference between MAGs and reference genomes. “SE”, “SA”, “PA”, “LM”, “LF”, “EF”, “EC” and “BS” is the abbreviation of bacterial species in Mock2 dataset. The assembly strategies have no significant influence on the genome size of MAGs. (C) The single-base error rate of MAGs compared with reference genomes. MAGs from MetaSPAdes assemblies have the highest single-base error rate of the four hybrid assembly strategies. The single base error rate of MAGs from long-read assemblies is significantly higher than that of hybrid assemblies. After being polished with long reads, the accuracy of MAGs has shown a slight improvement. Comparatively, the accuracy of MAGs is significantly improved after short-read and long-read polishing. (D) The number of breakpoints found in different MAGs. Long-read assembly strategies could produce more complete MAGs compared with hybrid assembly strategies. (E) The mean predicted gene length in different MAGs. Due to the high error rate, the mean gene length in MAGs from long-read assemblies is significantly shorter than in MAGs from hybrid assemblies. After short-read and long-read polishing, the mean predicted gene length is comparable to the mean gene length in MAGs from hybrid assemblies.


### Citation
Peng K, Gao Y, Li C, Wang Q, Yin Y, Hameed MF, Feil E, Chen S, Wang Z, Liu YX, Li R. Benchmarking of analysis tools and pipeline development for nanopore long-read metagenomics. Sci Bull (Beijing). 2025 Mar 20:S2095-9273(25)00310-X. doi: 10.1016/j.scib.2025.03.044.

Copyright 2023-2026 Kai Peng 008719@yzu.edu.cn
