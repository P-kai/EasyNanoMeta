[TOC]

# EasyNanoMeta:A pipeline for nanopore long-read metagenomic analysis

With the widespread application of nanopore sequencing in microbiology, there has been a significant increase in the generation of long-read sequencing data. This poses a challenge in effectively mining valuable biological information from the vast amount of nanopore sequencing data using bioinformatics tools. In this review, we present the advantages and pitfalls of long-read nanopore sequencing technologies in microbial genome and microbiome research, systematically summarize the bioinformatics analysis methods for nanopore long-read metagenomics and validate the practical implementation of the analysis workflow using various simulated and experimental metagenomic datasets. In addition, we conducted a comprehensive evaluation of the overall performance and computational requirement of various bioinformatic tools.

## Workflow overview

![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure1.jpg)

This comprehensive statement encompasses the various steps involved in the analysis process, along with the corresponding tools employed. The initial step in analyzing nanopore sequencing data involves data quality control, including basecalling, filtering and host sequences removal. Then, the the main steps of nanopore metagenomic data analysis were summarized into assembly-free strategy and assembly-based strategy. Assembly-free long-read metagenomic profiling can quickly estimate the abundance of species and functional genes employing few computing resources. In comparison, assembly-based nanopore metagenomic data analysis enable us to study more high-quality MAGs and their functions.

## How to use?  
1. Readme.md: Introduction  
2. install.sh: Dependencies installation. Using the script to install tools for nanopore metagenomic analysis.  
3. pipeline.sh: Command-line analysis for Linux. A guidance for analyzing nanopore metagenomic data.
4. Python scripts for data analysis: We provided many scripts for data processing.

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

``$ python co-located.py --help``  
``usage: co-located.py [-h] [--i I]``  
  
``options:``  
``  -h, --help  show this help message and exit``  
``  --i I``  
  
``$ python3 merge_tables.py -h``  
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
``
``Extracting_MAGs.py -h``  
``usage: test2.py [-h] [--flye FLYE] [--threads THREADS] --rawlong RAWLONG --rawshort1 RAWSHORT1 --rawshort2 RAWSHORT2``  
``                [--nextpolish NEXTPOLISH] [--semibin SEMIBIN] [--checkm2 CHECKM2] [--checkm2_database CHECKM2_DATABASE]``  
``                [--gtdbtk GTDBTK] [--gtdbtk_database GTDBTK_DATABASE]``  
``  
``Main script to run multiple bioinformatics tools.``  
``
``options:``  
``  -h, --help            show this help message and exit``  
``  --flye FLYE           Path to Flye executable (default: /Tools/software/Flye-2021-2.9/bin/flye).``  
``  --threads THREADS     Number of threads for MetaFlye.``  
``  --rawlong RAWLONG     Input raw long-read file.``  
``  --rawshort1 RAWSHORT1``  
``                        Input raw short-read1 file.``  
``  --rawshort2 RAWSHORT2``  
``                        Input raw short-read2 file.``  
``  --nextpolish NEXTPOLISH``  
``                        Path to NextPolish executable (default: /home/tools_pk/tools/NextPolish/nextPolish)``  
``  --semibin SEMIBIN     Path to SemiBin executable (default: /home/tools_pk/miniconda3/envs/SemiBin/bin/SemiBin).``  
``  --checkm2 CHECKM2     Path to the CheckM2 executable (default: /home/tools_pk/miniconda3/envs/checkm2/bin/checkm2).``  
``  --checkm2_database CHECKM2_DATABASE``  
``                        Path to the CheckM2 database (default:``  
``                        /home/tools_pk/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd).``  
``  --gtdbtk GTDBTK       Path to the GTDB-Tk executable (default: /home/tools_pk/miniconda3/envs/gtdbtk/bin/gtdbtk).``  
``  --gtdbtk_database GTDBTK_DATABASE``  
``                        Path to the GTDB-Tk database (default: /backup/database/gtdbtk/release214).``  

6. Rscripts_for_ploting: Examples for ploting using R.    

## Dependencies
**Miniconda3** (https://docs.anaconda.com/free/miniconda/index.html)  
**Minimap2** (https://github.com/lh3/minimap2)  
**samtools** (https://github.com/samtools/samtools)  
**bedtools** (https://github.com/arq5x/bedtools2)  
**fastp** (https://github.com/OpenGene/fastp)  
**Centrifuge** (http://www.ccb.jhu.edu/software/centrifuge)  
**Kraken2** (https://github.com/DerrickWood/kraken2)  
**Abricate** (https://github.com/tseemann/abricate)  
**MetaFlye** (https://github.com/fenderglass/Flye)  
**Canu** (https://github.com/marbl/canu)  
**wtdbg2** (https://github.com/ruanjue/wtdbg2)  
**NextDenove** (https://github.com/Nextomics/NextDenovo)  
**OPERA-MS** (https://github.com/CSB5/OPERA-MS)  
**MetaSPAdes** (https://github.com/ablab/spades)  
**MetaPlatanus** (https://github.com/rkajitani/MetaPlatanus)  
**Unicycler** (https://github.com/rrwick/Unicycler)  
**Racon** (https://github.com/isovic/racon)  
**NextPolish** (https://github.com/Nextomics/NextPolish)  
**Pilon** (https://github.com/broadinstitute/pilon)  
**bwa** (https://github.com/lh3/bwa)  
**Semibin** (https://github.com/BigDataBiology/SemiBin)  
**vamb** (https://github.com/RasmussenLab/vamb)  
**MetaWrap** (https://github.com/bxlab/metaWRAP)  
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
We used four hybrid assembly strategies and four long-read assembly strategies to perform metagenomic assembly of a Mock2 dataset. The long-read assemblies were polished with long-read data and combined with long-read and short-read data. The suffix “lp” indicates long-read polishing. The suffix “slp” indicates short-read and long-read polishing. After that, we extracted MAGs from different assemblies and compared MAGs with reference genomes. (A) Workflow for assessing the quality of MAGs. We compared the MAGs recovered using different analysis strategies with reference genomes to assess the quality of MAGs. (B) The genome size difference between MAGs and reference genomes. “SE”, “SA”, “PA”, “LM”, “LF”, “EF”, “EC” and “BS” is the abbreviation of bacterial species in Mock2 dataset (Table S4). The assembly strategies have no significant influence on the genome size of MAGs. (C) The single-base error rate of MAGs compared with reference genomes. MAGs from MetaSPAdes assemblies have the highest single-base error rate of the four hybrid assembly strategies. The single base error rate of MAGs from long-read assemblies is significantly higher than that of hybrid assemblies. After being polished with long reads, the accuracy of MAGs has shown a slight improvement. Comparatively, the accuracy of MAGs is significantly improved after short-read and long-read polishing. (D) The number of breakpoints found in different MAGs. Long-read assembly strategies could produce more complete MAGs compared with hybrid assembly strategies. (E) The mean predicted gene length in different MAGs. Due to the high error rate, the mean gene length in MAGs from long-read assemblies is significantly shorter than in MAGs from hybrid assemblies. After short-read and long-read polishing, the mean predicted gene length is comparable to the mean gene length in MAGs from hybrid assemblies.
