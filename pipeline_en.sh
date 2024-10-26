# EasyNanoMeta: Third-generation metagenomic analysis workflow 

    # Version: 1.01, 2023/7/20
    # Tested on CentOS 7.7

# 1. Acquisition of test data
## 1.1 Download raw data (mock dataset, activated sludge from sewage treatment plant, human gut metagenome, chicken gut metagenome)

### 1.1.1 Download mock dataset, reference: https://academic.oup.com/gigascience/article/8/5/giz043/5486468?login=true

    # Download third-generation metagenomes of mock dataset(including 8 bacteria and 2 yeasts, with a low proportion of yeast)
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/ERR003/3152/ERR3152364/ERR3152364.lite.1

    # Download second-generation metagenomes of mock dataset
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-16/ERR002/984/ERR2984773.sralite.1

    # Alternative download link of mock dataset(third-generation)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR3152364/ERR3152364
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR315/ERR3152364/Zymo-GridION-EVEN-BB-SN-flipflop.fq.gz

    # Alternative download link of mock dataset(second-generation)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR2984773/ERR2984773
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR298/ERR2984773/in732_1_R1.fastq.gz
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR298/ERR2984773/in732_1_R2.fastq.gz

### 1.1.2 Download human gut metagenome, reference: https://www.nature.com/articles/s41467-021-27917-x

    # Download third-generation human gut metagenome
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-zq-20/SRR017/17456/SRR17456359/SRR17456359.lite.1

    # Download second-generation human gut metagenome
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-zq-38/SRR017/17049/SRR17049035/SRR17049035.lite.1

    # Alternative download link of human gut metagenome(third-generation)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17456359/SRR17456359

    # Alternative download link of human gut metagenome(second-generation)
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17049035/SRR17049035

### 1.1.3 Download metagenome of activated sludge from sewage treatment plant, reference: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0663-0#Sec2

    # Download third-generation metagenome of activated sludge from sewage treatment plant
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-zq-22/SRR007/7627/SRR7627523/SRR7627523.lite.1

    # Download second-generation metagenome of activated sludge from sewage treatment plant
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR008/208/SRR8208348.sralite.1

    # Alternative download link of metagenome of activated sludge from sewage treatment plant(third-generation)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7627523/SRR7627523

    # Alternative download link of metagenome of activated sludge from sewage treatment plant(second-generation)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8208348/SRR8208348

## 1.2 Convert SRA data to fastq data

    # Use the variable i to specify the sample name
    i=sample_name

    # Decompress files with fastq-dump
    ~/tools/sratoolkit.3.0.6-centos_linux64/bin/fastq-dump --split-3 ${i}

    # Batch process samples, generate scripts first, then batch data conversion
    ls > list && for i in `cat list`; \
    do echo "~/tools/sratoolkit.3.0.6-centos_linux64/bin/fastq-dump --split-3 ${i}"; \
    done > fastq-dump.sh
    sh fastq-dump.sh

# 2. Quality control and filtering of test data

## 2.1 Quality control of third-generation raw data

### 2.1.1 Adapters removal of third-generation raw data
    # Adapter removal of nanopore data with porechop_abi
    # Activate software environment
    conda activate porechop_abi
    # Set the input sample name
    i=sample_name
    # Remove adapters with default parameters
    porechop_abi --ab_initio \
      -i ${i}.fastq \
      -o ${i}_output.fastq \
      -t 24

    # Batch process samples
    # Generate names for all samples
    ls *1.fastq && cut -f1 -d '.' > samples_name
    # Generate a script for batch processing
    for i in `cat samples_name`; \
    do echo "porechop_abi --ab_initio -1 ${i}.fastq -o ${i}_output.fastq -t 24"; \
    done > porechop_abi.sh
    # Run the script for batch processing
    sh porechop_abi.sh

### 2.1.2 Data distribution of quality control of third-generation raw data

    # Quality control of single sample data using NanoPlot
    # Activate the conda environment with NanoPlot installed
    conda activate nanoplot
    
    i=sample_name
    NanoPlot \
      -t 24 --N50 --huge -f svg --dpi 500 \
      --fastq ${i}.fastq \
      --plots kde dot \
      --title ${i} \
      -o ${i}

    # Quality control of multiple samples data using NanoPlot
    i1=sample_name1
    i2=sample_name2
    NanoComp --dpi 500 -t 24 -p prefix -f svg \
      --fastq ${i1}.fastq ${i2}.fastq \
      --names ${i1} ${i2} \
      -o NanoComp

### 2.1.3 Host removal of third-generation raw data, using minimap2, samtools and bedtools in combination

    # Use minimap2 and samtools for host removal of third-generation data, as an example of human gut microbiome removal of human genomic contamination
    # Download the human reference genome to the database folder
    cd ~/db
    wget https://ftp.ncbi.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/

    # Activate conda environment for host removal
    conda activate host_removal
    
    # Establish the minimap2 alignment index, with the human genome as the reference alignment sequence
    i=human_genome
    minimap2 -d ${i}.min ${i}.fasta

    # Use minimap2 for data comparison
    minimap2 -ax map-ont -t 24 ${i}.min ../raw.fasta -o minimap.sam

    # Extraction of sequences not matched to the host
    samtools view -bS -T -@24 ${i}.fasta -f 4 minimap.sam > unmaped_minimap.bam

    # Convert the bam file to a fastq file. First, sort the bam file, and then use the bamtofastq in bedtools to convert the bam file to the fastq file
    samtools sort -n unmaped_minimap.bam -o unmaped_sorted_minimap.bam
    bedtools bamtofastq -i unmaped_sorted_minimap.bam -fq fitted_raw.fastq

## 2.2 Quality control of second-generation raw data

### 2.2.1 Quality control of second-generation raw data
    
    # Use fastp to filter second-generation raw data
    i=sample_name
    fastp -i ${i}_raw_1.fq.gz \
      -I ${i}_raw_2.fq.gz \
      -o ${i}_clean_1.fq.gz \
      -O ${i}_clean_2.fq.gz \
      -t 24

    # Batch process samples using fastp
    # Generate names for all samples
    ls *1.fastq && cut -f1 -d '.' > samples_name
    # Generate a script for batch processing
    for i in `cat samples_name`; \
    do echo "fastp -i ${i}_1.fastq -I ${i}_2.fastq -o ${i}_clean_1.fastq -O ${i}_clean_2.fastq -t 24"; \
    done > fastp.sh
    # Run the script for batch processing
    sh porechop_abi.sh

# 3. Species and functional annotations of third-generation metagenomic raw data

## 3.1 Species annotations of third-generation metagenome

### 3.1.1 Use Centrifuge for species annotations of third-generation metagenome
    # Microbial species composition analysis using centrifuge with fastq data as input file
    i=sample_name
    ~/Pengkai/Tools/centrifuge/bin/centrifuge -p 24 \
      -x /path/to/Database/centrifuge_h+p+v_20200318/hpv \
      -q ${i}.fastq \
      --report-file ${i}_report \
      -S ${i}_result

    # Microbial species composition analysis using centrifuge with fasta data as input file
    i=sample_name
    ~/Pengkai/Tools/centrifuge/bin/centrifuge -p 24 \
      -x /path/to/Database/centrifuge_h+p+v_20200318/hpv \
      -f ${i}.fasta \
      --report-file ${i}_report \
      -S ${i}_result

    # Convert the centrifuge output to kraken2 results
    ~/Pengkai/Tools/centrifuge/bin/centrifuge-kreport \
      -x /ifs1/Database/centrifuge_h+p+v_20200318/hpv \
      ${i}_result > ${i}_kraken_report

### 3.1.2 Use kraken2 for species annotations of third-generation metagenome
    # Use kraken2 for species classification 
    i=sample_name
    kraken2 \
      --db /path/to/kraken2_db/k2_standard/ \
      --threads 24 \
      --report ${i}_kraken2_report \
      --output ${i}_kraken_result \
      ${i}.fastq

## 3.2 Annotations of antibiotic resistance genes in the third-generation metagenome

### 3.2.1 Use abricate to identify antibiotic resistance genes in the third-generation metagenomic raw data 
    # use abricate to identify antibiotic resistance genes in the third-generation metagenome
    i=sample_name
    abricate --db ncbi -t 24 ${i}.fasta > ${i}_ncbi_result

    # Use scripts to calculate the abundance of antibiotic resistance genes in metagenomes
    python abundance_calculate.py --i ${i}_ncbi_result --data_size metagenome_size > ${i}_abundance_result

### 3.2.2 Use abricate to identify insertion sequence in the third-generation metagenomic raw data 
    # Activate conda environment
    conda activate abricate

    # Download the latest insertion sequence database and build the local database, download address£ºhttps://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.fna
    # Enter the abricate database directory, and if it is installed in Conda or Miniconda, then in the Conda directory
    cd /path/to/abricate/db
    # Create an ISfinder folder and copy the fasta format database into the folder, rename it to sequences
    mkdir ISfinder
    cp /your/database/database.fasta ISfinder/sequences
    # Create a blast database named ISfinder
    makeblastdb -in sequences -title ISfinder -dbtype nucl -hash_index

    # Use abricate to identify insertion sequecnces in metagenome
    i=sample_name
    abricate --db ISfinder -t 24 ${i}.fasta > ${i}_ISfinder_result

### 3.2.3 Use abricate to identify virulence genes in the third-generation metagenomic raw data
    # Abricate provides a vfdb database. Update it
    abricate-get_db --db vfdb

    # Use abricate to identify insertion sequecnces in metagenome
    i=sample_name
    abricate --db vfdb -t 24 ${i}.fasta > ${i}_vfdb_result

# 4. Assembly of test data

## 4.1 Assembly of third-generation metagenome

### 4.1.1 Use flye for assembly of third-generation metagenome
    # Use the meta parameters of flye for assembly of third-generation metagenome
    i=sample_name
    # Software parameters
    ~/path/to/Flye-2.9.2/bin/flye \
      --meta \
      --nano-raw ${i}.fasta \
      --threads 24 \
      --out-dir ${i}_flye 

### 4.1.2 Use canu for assembly of third-generation metagenome
    # Use canu for assembly of long-read sequences
    i=sample_name
    /ifs1/User/yongxin/Pengkai/Tools/canu-2.2/bin/canu \
      -nanopore ${i}.fastq \
      -d canu \
      -p test \
      genomeSize=500m maxThreads=24
      useGrid=false # Solution to cluster operation failure

    # Use canu to correct errors of long-read data
    /ifs1/User/yongxin/Pengkai/Tools/canu-2.2/bin/canu \
      -correct \
      -nanopore half_output_reads_1.fastq \
      -d canu_correct \
      -p test \
      genomeSize=1g maxThreads=20 

### 4.1.3 Use wtdbg for assembly of third-generation metagenome
    # Software usage and parameters
    i=sample_name
    mkdir wtdbg2
    /ifs1/User/yongxin/Pengkai/Tools/wtdbg2/wtdbg2.pl -t 24 \
      -x ont \
      -g 200m \
      -o wtdbg2/CE4_1 \
      ${i}.fasta

### 4.1.4 Use NextDenovo for assembly of third-generation metagenome
    # To use the software, first prepare the input file input.fofn, and then edit the software running script run.cfg
    ls reads1.fasta > input.fofn
    vim run.cfg

        ###
        [General]
        job_type = local
        job_prefix = nextDenovo
        task = all # 'all', 'correct', 'assemble'
        rewrite = yes # yes/no
        deltmp = yes
        rerun = 3
        parallel_jobs = 2
        input_type = raw
        read_type = clr
        input_fofn = ./input.fofn
        workdir = ./01_rundir

        [correct_option]
        read_cutoff = 300bp # Set sequence filtering length
        genome_size = 300008161 # Set genome assembly size
        pa_correction = 2
        sort_options = -m 1g -t 2
        minimap2_options_raw =  -t 8
        correction_options = -p 15

        [assemble_option]
        minimap2_options_cns =  -t 8
        nextgraph_options = -a 1
        ###

    # Run assembly script
    /ifs1/User/yongxin/Pengkai/Tools/NextDenovo/nextDenovo run.cfg 

## 4.2 Hybrid assembly of second-generation and third-generation metagenomes

### 4.2.1 Use OPERA-MS for hybrid assembly of second-generation and third-generation metagenomic data

    # Activate operams environment, Use test data for software testing
    conda activate operams
    perl ../OPERA-MS.pl \
        --contig-file contigs.fasta \
        --short-read1 R1.fastq.gz \
        --short-read2 R2.fastq.gz \
        --long-read long_read.fastq \
        --no-ref-clustering \
        --num-processors 24 \
        --out-dir RESULTS       

    # Use OPERA-MS for hybrid assembly of second-generation and third-generation metagenomes with spades as assembly mode of second-generation sequences
    conda activate operams
    perl ../OPERA-MS.pl \
        --short-read1 R1.fastq.gz \
        --short-read2 R2.fastq.gz \
        --long-read long_read.fastq \
        --no-ref-clustering \
        --num-processors 24 \
        --out-dir RESULTS

### 4.2.2 Use metaSPAdes for hybrid assembly of second-generation and third-generation metagenomic data

    # Use spades for hybrid assembly of second-generation and third-generation metagenomes
    ~/Pengkai/Tools/SPAdes-3.15.5-Linux/bin/spades.py --meta \
      -1 SRR17049035.lite.1_1.fastq -2 SRR17049035.lite.1_2.fastq \
      --nanopore SRR17456359.lite.1.fastq \
      -t 24 -o metaSPAdes_human_C29 \
      --phred-offset 33

### 4.2.3 Use unicycler for hybrid assembly of second-generation and third-generation metagenomic data

    # Use unicycler for hybrid assembly of second-generation and third-generation metagenomes
    i=sample_name

    unicycler -1 ${i}_1.fastq -2 ${i}_2.fastq \
    -l ${i}.fastq \
    -o ${i}_unicycler \
    -t 24

### 4.2.4 Use metaplatanus for hybrid assembly of second-generation and third-generation metagenomic data
   
    # Use metaplatanus for hybrid assembly of second-generation and third-generation metagenomes
    i=sample_name

    metaplatanus -IP2 ${i}_1.fastq ${i}_2.fastq \
    -ont ${i}.fastq \
    -t 24 \
    -o ${i}_metaplatanus

# 5. Calibrate third-generation metagenome assembly result

## 5.1 Calibrate assembly result of lonely third-generation metagenome 

### 5.1.1 Use minimap2 and racon to calibrate third-generation metagenome assembly result

    # Index the assembled genomes
    i=sample_name
    minimap2 -d ref_1.mmi ${i}.fa

    # Align sequences with minimap2. raw.fasta is unassembled raw sequencing data, also available in fastq format
    minimap2 -ax map-ont -t 8  ref_1.mmi raw.fasta > alin_1.sam
    
    # Calibrate assembled genomes with racon
    racon -t 8 raw.fasta alin_1.sam ${i}.fasta > ${i}_polished1.fa  

### 5.1.2 Use NextPolish to calibrate third-generation metagenome assembly result, based on short-read high-accuracy second-generation metagenomic data

    # Calibrate third-generation metagenome assembly result with second-generation metagenomic data
    # Generate file locations of short-read metagenomic data for calibration of third-generation metagenome assembly result
    ls reads1.fq reads2.fa.gz > sgs.fofn 
    
    # Edit executable file of assembly result calibration and configure the software execution parameters
    vim run.cfg  

        ###
        [General]
        job_type = local
        job_prefix = nextPolish
        task = best
        rewrite = yes
        rerun = 3
        parallel_jobs = 6
        multithread_jobs = 5
        genome = ./raw.genome.fasta # Assembly result file
        genome_size = auto
        workdir = ./short-reads-polish
        polish_options = -p 8

        [sgs_option]
        sgs_fofn = ./sgs.fofn
        sgs_options = -max_depth 100 -bwa
        ### Edit the above information in the run.cfg file

    # Perform calibration procedures for assembly result
    path/nextPolish run.cfg

### 5.1.3 Use NextPolish to calibrate third-generation metagenome assembly result, based on third-generation metagenomic raw data

    # Calibrate third-generation metagenome assembly result with third-generation metagenomic raw data
    # Generate file locations of long-read metagenomic data for calibration of third-generation metagenome assembly result
    ls raw_nanopore.fq > lgs.fofn 
    
    # Edit executable file of assembly result calibration and configure the software execution parameters
    vim run.cfg  

        ###
        [General]
        job_type = local
        job_prefix = nextPolish
        task = best
        rewrite = yes
        rerun = 3
        parallel_jobs = 6
        multithread_jobs = 5
        genome = ./raw.genome.fasta # Assembly result file
        genome_size = auto
        workdir = ./01_rundir
        polish_options = -p 8

        [lgs_option]
        lgs_fofn = ./lgs.fofn
        lgs_options = -min_read_len 1k -max_depth 100
        lgs_minimap2_options = -x map-ont
        ### Edit the above information in the run.cfg file

    # Perform calibration procedures for assembly result
    path/nextPolish run.cfg

### 5.1.4 Use NextPolish to calibrate third-generation metagenome assembly result, based on second-generation and third-generation metagenomic raw data

    # Calibrate third-generation metagenome assembly result with long-read and short-read metagenomic raw data
    # Generate file locations of long-read and short-read metagenomic data for calibration of third-generation metagenome assembly result
    ls reads1.fq reads2.fa.gz > sgs.fofn
    ls raw_nanopore.gz > lgs.fofn

    # Edit executable file of assembly result calibration and configure the software execution parameters
    vim run.cfg  

        ###
        [General]
        job_type = local
        job_prefix = nextPolish
        task = best
        rewrite = yes
        rerun = 3
        parallel_jobs = 6
        multithread_jobs = 5
        genome = ./raw.genome.fasta  # Assembly result file
        genome_size = auto
        workdir = ./short-and-long-reads-polish
        polish_options = -p 8

        [sgs_option]
        sgs_fofn = ./sgs.fofn
        sgs_options = -max_depth 100 -bwa

        [lgs_option]
        lgs_fofn = ./lgs.fofn
        lgs_options = -min_read_len 1k -max_depth 100
        lgs_minimap2_options = -x map-ont
        ### Edit the above information in the run.cfg file

    # Perform calibration procedures for assembly result
    path/nextPolish run.cfg

### 5.1.5 Use Pilon to calibrate third-generation metagenome assembly result, based on second-generation metagenomic raw data

    # Start calibration of assembly result
    # Index assembly result using bwa
    i=sample_name
    bwa index ${i}.fasta
    # Map second-generation short-read metagenomic data to assembly results using bwa
    bwa mem -t 6 ${i}.fasta raw_read1.fq raw_read2.fq > ${i}.sam
    # Convert sam file to bam file
    samtools view -@ 10 -bS -F 4 ${i}.sam > ${i}.bam
    # Sort bam file
    samtools sort -@ 10 ${i}.bam > ${i}.sorted.bam
    # Index sorted bam file
    samtools index ${i}.sorted.bam
    # Use pilon for assembly result calibration
    java -Xmx16G -jar path/to/pilon-1.22.jar --genome ${i}.fasta --frags aln-pe.sorted.bam --output ${i}_polished

    # Carefully, calibration of assembly result using pilon generally requires multiple rounds of calibration, 3-4 rounds recommended.

# 6. Binning and binning purify and reassembly analysis of third-generation metagenome assembly result

## 6.1 Bin assembly result of lonely third-generation metagenome 

### 6.1.1 Use SemiBin for binning of third-generation metagenome assembly result

    # Software usage
    # Index assembled genome
    i=sample_name
    minimap2 -d catalogue.mmi ${i}.fasta
    # Align index file with raw sequencing data to get bam file, raw_data.fq.gz is raw sequencing data
    minimap2 -t 8 -N 5 -ax map-ont catalogue.mmi --split-prefix mmsplit ../raw_data.fq.gz | samtools view -F 3584 -b --threads 8 > ${i}.bam
    # Sort bam file
    samtools sort -@ 10 ${i}.bam > ${i}.sorted.bam 
    # Use SemiBin to run binning
    SemiBin single_easy_bin -i ${i}.fasta  --sequencing-type long_read -b ${i}.sorted.bam -o bin_output --environment global 

## 6.2 Bin hybrid assembly result of second-generation and third-generation metagenomes

### 6.2.1 Use metawrap for binning of second-generation and third-generation metagenomes hybrid assembly result
    # Use binning mode of metawrap for binning of second-generation and third-generation metagenomes hybrid assembly result
    # Run binning with second-generation sequencing data as raw data, format must be *_1.fastq and *_2.fastq
    i=sample_name
    metawrap binning --metabat2 --maxbin2 --concoct -t 48 --run-checkm -a ${i}.fa -o bin ${i}_clean_1.fastq ${i}_clean_2.fastq

 ## 6.3 Binning purify and quantification of third-generation metagenome assembly result
    # Software£ºmetawrap
    # Use metawrap for Binning purification, with parameters of integrity greater than 80% and contamination less than 10%
    metawrap bin_refinement \
      -A maxbin2_bins/ \
      -B metabat2_bins/ \
      -C concoct_bins/ \
      -o out_dir \
      -t 12 -c 80 

    # Binning quantification based on second-generation metagenomic raw data
    metawrap quant_bins -t 24 \
      -b metawrap_80_10_bins \
      -o bin_quant \
      -a ${i}.fa ${i}.fa_clean_1.fastq ${i}.fa_clean_2.fastq 

    # Use vamb for binning
    vamb --outdir canu_vamb --fasta canu.contigs.fasta --bamfiles canu.sorted.bam --minfasta 200000

# 7. Species and functional annotations of binning result

## 7.1 Species annotations of binning result MAG

### 7.1.1 Use GTDB-Tk for species annotations of MAG
    # MAG genomic classifications and annotations software GTDB-Tk
    # Use conda to create new environment and install GTDB-Tk
    conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1 
    # Activate GTDB-Tk environment
    conda activate gtdbtk-2.1.1 
    # Directly use the script that comes with the software to configure the database, the database configuration location is ~/miniconda3/envs/gtdbtk-2.1.1/share/gtdbtk-2.1.1/db/gtdbtk_r207_v2_data.tar.gz
    download-db.sh

    # Option to manually configure the database is also available
    # Set the database path
    mkdir -p ~/db/gtdb & cd ~/db/gtdb
    # Download and decompress database
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz  
    tar -zxvf gtdbtk_data.tar.gz
    # Set the database location, pay attention to modify the software installation location
    locate gtdbtk.sh # Locate configured file
    # Modify the path after PATH= to the directory where the database was extracted, like /home/meta/db/gtdb/release95/
    vim /conda/envs/gtdbtk/etc/conda/activate.d/gtdbtk.sh

    # Download the appropriate database from the official website, here download the 214 version of the database
    wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz 
    # Decompress database
    tar -zcvf gtdbtk_r214_data.tar.gz
    # Configure the downloaded database to the software environment variable
    export GTDBTK_DATA_PATH=/ifs1/User/yongxin/db/gtdb 

    # Determine whether third-party programs and databases required by the software are working properly
    gtdbtk check_install 
    # Run binning species classifications and construct evolutionary tree
    gtdbtk classify_wf --genome_dir bins/ --extension fa  --skip_ani_screen --out_dir gtdbtk 
    gtdbtk convert_to_itol --input some_tree.tree --output itol.tree # Convert evolutionary tree file to itol format

## 7.2 Functional annotations of binning result MAG
    # Use prokka for functional annotations of binning result MAG
    i=sample_name
    prokka ${i}.fa --prefix ${i} --outdir ~/prokka/${i}

    # Use abricate to annotate antibiotic resistance genes in MAG
    abricate --db ncbi --minid 80 --mincov 80 -t 12 *fa > ncbi_annotation_results
