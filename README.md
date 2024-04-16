[TOC]

# EasyNanoMeta:A pipeline for nanopore long-read metagenomic analysis

With the widespread application of nanopore sequencing in microbiology, there has been a significant increase in the generation of long-read sequencing data. This poses a challenge in effectively mining valuable biological information from the vast amount of nanopore sequencing data using bioinformatics tools. In this review, we present the advantages and pitfalls of long-read nanopore sequencing technologies in microbial genome and microbiome research, systematically summarize the bioinformatics analysis methods for nanopore long-read metagenomics and validate the practical implementation of the analysis workflow using various simulated and experimental metagenomic datasets. In addition, we conducted a comprehensive evaluation of the overall performance and computational requirement of various bioinformatic tools.

## Workflow overview

![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure2.png)

The main analysis pipeline for nanopore microbial sequences. This comprehensive statement encompasses the various steps involved in the analysis process, along with the corresponding tools employed. The initial step in analyzing nanopore sequencing data involves data quality control, including basecalling, filtering and host sequences removal. Presently, several tools can be utilized for the analysis of full-length 16S sequencing data obtained through nanopore sequencing. Emu is the most reliable and user-friendly option in this regard. For bacterial genome analysis, accurate and complete bacterial genomes can be generated via two distinct strategies. The first approach entails performing long-read assembly and subsequently refining the assembly using high-accuracy short-read data. The second approach involves hybrid assembly, which combines both long-read and short-read data. Nanopore metagenomic data analysis is more complex than 16S data and bacterial genome data. Assembly-free long-read metagenomic profiling can quickly estimate the abundance of species and functional genes employing few computing resources. In comparison, assembly-based nanopore metagenomic data analysis enable us to study more high-quality MAGs and their functions.

## Pipeline manual and file description  
1. Readme.md: Introduction  
2. install.sh: Dependencies installation. Using the script to install tools for nanopore metagenomic analysis.  
3. pipeline.sh: Command-line analysis for Linux. A guidance for analyzing nanopore metagenomic data.
4. Python scripts for data analysis: We provided many scripts for data processing.
5. Rscripts_for_ploting: Examples for ploting using R.    

## A practical guidance for assembly-free long-read metagenomic data analysis

![image](https://github.com/P-kai/EasyNanoMeta/blob/main/Figures/figure3.png)

We used a set of datasets from WWTP to perform assembly-free analysis. After obtaining the taxonomic and resistance abundance tables, we can perform many statistical and visualization analyses, including composition, diversity, and correlation. (A) The workflow for assembly-free nanopore metagenomic data analysis. (B) The microbial compositions at the phylum level. The most abundant bacteria in AS and SIF were Proteobacteria. (C) The distribution of the top 20 abundant ARGs in different samples. ARGs in SIF samples are more abundant than in AS samples. (D) Alpha diversity analysis for microbial and ARG compositions. Alpha diversity showed that AS group has more abundant microbes but less abundant ARGs. (E) Beta diversity analysis for ARG compositions. PCoA based on Bray-Curtis showed that there is a significant difference in ARG compositions between AS and SIF. (F) The correlation analysis of microbial and ARG compositions. In the heat map, dark blue color blocks indicate positive correlations between different ARGs based on the Pearson Correlation Coefficient. The attachment of ARGs and the bacterial phylum indicate the correlation. (G) The correlation analysis of different ARGs based on their genetic context. The co-occurrence network diagram shows the correlation of different ARGs based on their genetic location. The nodes represent different ARGs. Connections between nodes indicate that they are related. The edge widths represent the correlation coefficient.

## Dependencies
