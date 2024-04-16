# 三代宏基因组分析流程EasyNanoMeta

    # 版本: 1.01, 2023/7/20
    # 测试环境为CentOS 7.7

# 一、测试数据获取
## 1.1 下载原始数据（mock数据集，污水处理厂活性污泥，人肠道宏基因组，鸡肠道宏基因组）

### 1.1.1 mock数据集下载，参考文献：https://academic.oup.com/gigascience/article/8/5/giz043/5486468?login=true

    # 下载mock数据集(包括8个细菌，2个酵母，酵母占比较低)的三代数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/ERR003/3152/ERR3152364/ERR3152364.lite.1

    # 下载mock数据集的二代数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-16/ERR002/984/ERR2984773.sralite.1

    # mock数据备用下载链接（三代）
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR3152364/ERR3152364
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR315/ERR3152364/Zymo-GridION-EVEN-BB-SN-flipflop.fq.gz

    # mock数据备用下载链接（二代）
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR2984773/ERR2984773
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR298/ERR2984773/in732_1_R1.fastq.gz
    wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR298/ERR2984773/in732_1_R2.fastq.gz

### 1.1.2 人肠道宏基因组数据下载，参考文献：https://www.nature.com/articles/s41467-021-27917-x

    # 下载人肠道三代宏基因组数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-zq-20/SRR017/17456/SRR17456359/SRR17456359.lite.1

    # 下载人肠道二代宏基因组数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-zq-38/SRR017/17049/SRR17049035/SRR17049035.lite.1

    # 人肠道宏基因组数据备用连接(三代)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17456359/SRR17456359

    # 人肠道宏基因组数据备用连接(二代)
    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR17049035/SRR17049035

### 1.1.3 污水处理厂活性污泥宏基因组数据，参考文献：https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0663-0#Sec2

    # 下载污水处理厂三代宏基因组数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-zq-22/SRR007/7627/SRR7627523/SRR7627523.lite.1

    # 下载污水处理厂二代宏基因组数据
    wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR008/208/SRR8208348.sralite.1

    # 污水处理厂活性污泥宏基因组数据备用连接(三代)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7627523/SRR7627523

    # 污水处理厂活性污泥宏基因组数据备用连接(二代)
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8208348/SRR8208348

## 1.2 SRA数据转换fastq数据

    # 使用变量i指定样本名称
    i=sample_name

    # 使用fastq-dump解压缩文件
    ~/tools/sratoolkit.3.0.6-centos_linux64/bin/fastq-dump --split-3 ${i}

    # 批量处理样本，先生成脚本，再进行批量数据转换
    ls > list && for i in `cat list`; \
    do echo "~/tools/sratoolkit.3.0.6-centos_linux64/bin/fastq-dump --split-3 ${i}"; \
    done > fastq-dump.sh
    sh fastq-dump.sh

# 二、对测试数据进行质控，过滤

## 2.1 三代原始数据质控

### 2.1.1 三代原始数据去接头
    # 使用porechop_abi进行nanopore数据接头去除
    # 激活软件所在环境
    conda activate porechop_abi
    # 设置输入样本名称
    i=sample_name
    # 使用默认参数去除接头
    porechop_abi --ab_initio \
      -i ${i}.fastq \
      -o ${i}_output.fastq \
      -t 24

    # 批量处理样本
    # 生成所有样本的名称
    ls *1.fastq && cut -f1 -d '.' > samples_name
    # 生成批量运行的脚本
    for i in `cat samples_name`; \
    do echo "porechop_abi --ab_initio -1 ${i}.fastq -o ${i}_output.fastq -t 24"; \
    done > porechop_abi.sh
    # 运行批量执行脚本
    sh porechop_abi.sh

### 2.1.2 三代原始数据质量控制数据分布

    # 使用NanoPlot进行单样本的数据质控
    # 激活进入安装NanoPlot的conda环境
    conda activate nanoplot
    
    i=sample_name
    NanoPlot \
      -t 24 --N50 --huge -f svg --dpi 500 \
      --fastq ${i}.fastq \
      --plots kde dot \
      --title ${i} \
      -o ${i}

    # 使用NanoComp进行多样本数据质控
    i1=sample_name1
    i2=sample_name2
    NanoComp --dpi 500 -t 24 -p prefix -f svg \
      --fastq ${i1}.fastq ${i2}.fastq \
      --names ${i1} ${i2} \
      -o NanoComp

### 2.1.3 三代原始数据去宿主，minimap2，samtools，bedtools连用

    #使用minimap2与samtools进行三代数据去宿主，以人肠道微生物组为例，进行人基因组污染去除
    #下载人的参考基因组到数据库文件夹下
    cd ~/db
    wget https://ftp.ncbi.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/

    # 进入去宿主的conda环境
    conda activate host_removal
    
    # 建立minimap2比对索引，human基因组为参考对齐序列
    i=human_genome
    minimap2 -d ${i}.min ${i}.fasta

    # 使用minimap2进行数据比对
    minimap2 -ax map-ont -t 24 ${i}.min ../raw.fasta -o minimap.sam

    # 提取未匹配到宿主的序列
    samtools view -bS -T -@24 ${i}.fasta -f 4 minimap.sam > unmaped_minimap.bam

    # 将bam文件转换为fastq文件，首先对bam文件进行排序，然后使用bedtools中的bamtofastq进行bam文件到fastq文件的转换
    samtools sort -n unmaped_minimap.bam -o unmaped_sorted_minimap.bam
    bedtools bamtofastq -i unmaped_sorted_minimap.bam -fq fitted_raw.fastq

## 2.2 二代原始数据质控

### 2.2.1 二代原始数据质量控制数据分布
    
    #使用fastp过滤二代原始数据
    i=sample_name
    fastp -i ${i}_raw_1.fq.gz \
      -I ${i}_raw_2.fq.gz \
      -o ${i}_clean_1.fq.gz \
      -O ${i}_clean_2.fq.gz \
      -t 24

# 三、三代宏基因组原始数据物种及功能注释

## 3.1 三代宏基因组数据物种注释

### 3.1.1 使用centrifuge进行三代宏基因组物种注释
    #使用centrifuge进行微生物物种组成分析，输入文件为fastq数据
    i=sample_name
    ~/Pengkai/Tools/centrifuge/bin/centrifuge -p 24 \
      -x /path/to/Database/centrifuge_h+p+v_20200318/hpv \
      -q ${i}.fastq \
      --report-file ${i}_report \
      -S ${i}_result

    #使用centrifuge进行微生物物种组成分析，输入文件为fasta数据
    i=sample_name
    ~/Pengkai/Tools/centrifuge/bin/centrifuge -p 24 \
      -x /path/to/Database/centrifuge_h+p+v_20200318/hpv \
      -f ${i}.fasta \
      --report-file ${i}_report \
      -S ${i}_result

    #将上一步centrifuge输出结果转换为kraken2结果
    ~/Pengkai/Tools/centrifuge/bin/centrifuge-kreport \
      -x /ifs1/Database/centrifuge_h+p+v_20200318/hpv \
      ${i}_result > ${i}_kraken_report

### 3.1.2 使用kraken2进行三代宏基因组物种注释
    #使用kraken2进行物种分类
    i=sample_name
    kraken2 \
      --db /path/to/kraken2_db/k2_standard/ \
      --threads 24 \
      --report ${i}_kraken2_report \
      --output ${i}_kraken_result \
      ${i}.fastq

## 3.2 三代宏基因组数据抗性基因注释

### 3.2.1 使用abricate进行三代宏基因组原始数据抗性基因识别
    #使用abricate识别宏基因组序列中的抗性基因
    i=sample_name
    abricate --db ncbi -t 24 ${i}.fasta > ${i}_ncbi_result

    #使用脚本进行宏基因组中耐药基因丰度计算
    python abundance_calculate.py --i ${i}_ncbi_result --data_size metagenome_size > ${i}_abundance_result

### 3.2.2 使用abricate进行三代宏基因组原始数据插入序列识别
    #激活conda环境
    conda activate abricate

    #下载最新的插入序列数据库并构建库，下载地址：https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.fna
    #进入abricate数据库目录，若为conda或miniconda安装，则在conda目录下
    cd /path/to/abricate/db
    #创建ISfinder文件夹，并将fasta格式的数据库拷贝到文件夹内，重命名为sequences
    mkdir ISfinder
    cp /your/database/database.fasta ISfinder/sequences
    #创建数据库，名称为 ISfinder
    makeblastdb -in sequences -title ISfinder -dbtype nucl -hash_index

    #使用abricate识别宏基因组序列中的插入序列
    i=sample_name
    abricate --db ISfinder -t 24 ${i}.fasta > ${i}_ISfinder_result

### 3.2.3 使用abricate进行三代宏基因组原始数据毒力基因的识别
    #软件提供vfdb毒力基因数据库，更新该数据库
    abricate-get_db --db vfdb

    #使用abricate识别宏基因组序列中的毒力基因
    i=sample_name
    abricate --db vfdb -t 24 ${i}.fasta > ${i}_vfdb_result

# 四、测试数据组装

## 4.1 三代宏基因组组装

### 4.1.1 使用flye进行三代宏基因组组装
    #使用flye的meta参数进行三代宏基因组组装
    i=sample_name
    #软件参数
    ~/path/to/Flye-2.9.2/bin/flye \
      --meta \
      --nano-raw ${i}.fasta \
      --threads 24 \
      --out-dir ${i}_flye 

### 4.1.2 使用canu进行三代宏基因组组装
    #使用canu进行长读数据组装
    i=sample_name
    /ifs1/User/yongxin/Pengkai/Tools/canu-2.2/bin/canu \
      -nanopore ${i}.fastq \
      -d canu \
      -p test \
      genomeSize=500m maxThreads=24
      useGrid=false #集群运行失败的解决办法

    #使用canu进行长读数据纠错
    /ifs1/User/yongxin/Pengkai/Tools/canu-2.2/bin/canu \
      -correct \
      -nanopore half_output_reads_1.fastq \
      -d canu_correct \
      -p test \
      genomeSize=1g maxThreads=20 

### 4.1.3 使用wtdbg2进行三代宏基因组组装
    #软件使用及参数
    i=sample_name
    mkdir wtdbg2
    /ifs1/User/yongxin/Pengkai/Tools/wtdbg2/wtdbg2.pl -t 24 \
      -x ont \
      -g 200m \
      -o wtdbg2/CE4_1 \
      ${i}.fasta

### 4.1.4 使用NextDenovo进行三代宏基因组组装
    #软件使用，首先准备输入文件input.fofn，再编辑软件运行脚本run.cfg
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
        read_cutoff = 300bp #设置序列过滤长度
        genome_size = 300008161 #设置基因组组装大小
        pa_correction = 2
        sort_options = -m 1g -t 2
        minimap2_options_raw =  -t 8
        correction_options = -p 15

        [assemble_option]
        minimap2_options_cns =  -t 8
        nextgraph_options = -a 1
        ###

    #运行组装脚本
    /ifs1/User/yongxin/Pengkai/Tools/NextDenovo/nextDenovo run.cfg 

## 4.2 二三代宏基因组混合组装

### 4.2.1 使用OPERA-MS进行二三代宏基因组数据组装

    #激活operams环境，使用测试数据进行软件测试
    conda activate operams
    perl ../OPERA-MS.pl \
        --contig-file contigs.fasta \
        --short-read1 R1.fastq.gz \
        --short-read2 R2.fastq.gz \
        --long-read long_read.fastq \
        --no-ref-clustering \
        --num-processors 24 \
        --out-dir RESULTS       

    #使用OPERA-MS进行二三代混合组装，二代组装模式选用spades
    conda activate operams
    perl ../OPERA-MS.pl \
        --short-read1 R1.fastq.gz \
        --short-read2 R2.fastq.gz \
        --long-read long_read.fastq \
        --no-ref-clustering \
        --num-processors 24 \
        --out-dir RESULTS

### 4.2.2 使用metaSPAdes进行二三代宏基因组数据组装

    #使用spades进行二三代混合组装
    ~/Pengkai/Tools/SPAdes-3.15.5-Linux/bin/spades.py --meta \
      -1 SRR17049035.lite.1_1.fastq -2 SRR17049035.lite.1_2.fastq \
      --nanopore SRR17456359.lite.1.fastq \
      -t 24 -o metaSPAdes_human_C29 \
      --phred-offset 33

### 4.2.3 使用unicycler进行二三代宏基因组数据组装

    #使用unicycler进行二三代混合组装
    i=sample_name

    unicycler -1 ${i}_1.fastq -2 ${i}_2.fastq \
    -l ${i}.fastq \
    -o ${i}_unicycler \
    -t 24

### 4.2.4 使用metaplatanus进行二三代宏基因组数据组装
   
    #使用metaplatanus进行二三代混合组装
    i=sample_name

    metaplatanus -IP2 ${i}_1.fastq ${i}_2.fastq \
    -ont ${i}.fastq \
    -t 24 \
    -o ${i}_metaplatanus

# 五、三代组装结果校准

## 5.1 纯三代宏基因组组装结果校准

### 5.1.1 使用minimap2与racon进行三代宏基因组数据组装结果校准

    #对组装基因组建立索引
    i=sample_name
    minimap2 -d ref_1.mmi ${i}.fa

    #使用minimap2对齐序列,raw.fasta为测序未组装原始数据，也可为fastq格式
    minimap2 -ax map-ont -t 8  ref_1.mmi raw.fasta > alin_1.sam
    
    #使用racon进行基因组校准
    racon -t 8 raw.fasta alin_1.sam ${i}.fasta > ${i}_polished1.fa  

### 5.1.2 以二代数据短读高准确数据为基础，使用NextPolish进行三代宏基因组数据组装结果校准

    #使用短读二代数据进行三代组装结果校准
    #生成短读数据文件位置，用于校准组装结果
    ls reads1.fq reads2.fa.gz > sgs.fofn 
    
    #编辑组装结果校准的可执行文件，在该文件中配置软件执行参数
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
        genome = ./raw.genome.fasta #组装结果文件
        genome_size = auto
        workdir = ./short-reads-polish
        polish_options = -p 8

        [sgs_option]
        sgs_fofn = ./sgs.fofn
        sgs_options = -max_depth 100 -bwa
        ### 在run.cfg文件中编辑以上信息

    #执行组装结果校准程序
    path/nextPolish run.cfg

### 5.1.3 以三代原始数据为基础，使用NextPolish进行三代宏基因组数据组装结果校准

    #使用三代原始数据进行三代组装结果校准
    #生成长读数据文件位置，用于校准组装结果
    ls raw_nanopore.fq > lgs.fofn 
    
    #编辑组装结果校准的可执行文件，在该文件中配置软件执行参数
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
        genome = ./raw.genome.fasta #组装结果文件
        genome_size = auto
        workdir = ./01_rundir
        polish_options = -p 8

        [lgs_option]
        lgs_fofn = ./lgs.fofn
        lgs_options = -min_read_len 1k -max_depth 100
        lgs_minimap2_options = -x map-ont
        ### 在run.cfg文件中编辑以上信息

    #执行组装结果校准程序
    path/nextPolish run.cfg

### 5.1.4 以二代、三代原始数据为基础，使用NextPolish进行三代宏基因组数据组装结果校准

    #使用长读、短读数据进行组装基因组校准，首先生成短读、长读数据文件位置
    ls reads1.fq reads2.fa.gz > sgs.fofn
    ls reads1.fq reads2.fa.gz > lgs.fofn

    #编辑组装结果校准执行文件，配置软件执行参数
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
        genome = ./raw.genome.fasta  #组装结果文件
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
        ### 在run.cfg文件中编辑以上信息

    #执行组装结果校准程序       
    path/nextPolish run.cfg

### 5.1.5 以二代原始数据为基础，使用Pilon进行三代宏基因组数据组装结果校准

    #开始进行组装结果校准
    #使用bwa对组装结果建立索引
    i=sample_name
    bwa index ${i}.fasta
    #使用bwa将二代短读数据映射到组装结果上
    bwa mem -t 6 ${i}.fasta raw_read1.fq raw_read2.fq > ${i}.sam
    #将sam文件转换为bam文件
    samtools view -@ 10 -bS -F 4 ${i}.sam > ${i}.bam
    #对bam文件进行排序
    samtools sort -@ 10 ${i}.bam > ${i}.sorted.bam
    #索引排序后的bam文件
    samtools index ${i}.sorted.bam
    #使用下载的pilon进行组装结果校准
    java -Xmx16G -jar path/to/pilon-1.22.jar --genome ${i}.fasta --frags aln-pe.sorted.bam --output ${i}_polished

    #注意，使用pilon校准组装结果，一般进行多轮校准，推荐3-4轮

# 六、三代宏基因组组装结果分箱及分箱提纯、重组装分析

## 6.1 纯三代宏基因组组装结果分箱

### 6.1.1 使用SemiBin进行三代宏基因组组装结果分箱

    #软件使用
    #首先对组装基因组进行索引创建，获取排序的bam文件
    i=sample_name
    minimap2 -d catalogue.mmi ${i}.fasta
    #比对获取bam文件,raw_data.fq.gz为测序的原始数据
    minimap2 -t 8 -N 5 -ax map-ont catalogue.mmi --split-prefix mmsplit ../raw_data.fq.gz | samtools view -F 3584 -b --threads 8 > ${i}.bam
    #对bam文件进行排序
    samtools sort -@ 10 ${i}.bam > ${i}.sorted.bam 
    #使用SemiBin运行bin
    SemiBin single_easy_bin -i ${i}.fasta  --sequencing-type long_read -b ${i}.sorted.bam -o bin_output --environment global 

## 6.2 二三代宏基因组混合组装结果分箱

### 6.2.1 使用metawrap进行二三代宏基因组组装结果分箱
    #使用metawrap里面的分箱模式进行二三代宏基因组组装数据分箱
    #运行分箱，原始数据为二代数据，格式必须为*_1.fastq；*_2.fastq
    i=sample_name
    metawrap binning --metabat2 --maxbin2 --concoct -t 48 --run-checkm -a ${i}.fa -o bin ${i}_clean_1.fastq ${i}_clean_2.fastq

 ## 6.3 宏基因组组装结果分箱提纯及定量
    #软件：metawrap
    #使用metawrap进行分箱提纯，参数为完整度大于80%，污染小于10%
    metawrap bin_refinement \
      -A maxbin2_bins/ \
      -B metabat2_bins/ \
      -C concoct_bins/ \
      -o out_dir \
      -t 12 -c 80 

    #分箱定量，此处为根据二代宏基因组原始数据进行定量
    metawrap quant_bins -t 24 \
      -b metawrap_80_10_bins \
      -o bin_quant \
      -a ${i}.fa ${i}.fa_clean_1.fastq ${i}.fa_clean_2.fastq 


#安装vamb
pip install vamb -i https://pypi.tuna.tsinghua.edu.cn/simple
/ifs1/User/yongxin/.local/bin #软件路径
export PATH=$PATH:/ifs1/User/yongxin/.local/bin/ #临时添加软件路径到环境变量
vamb --outdir canu_vamb --fasta canu.contigs.fasta --bamfiles canu.sorted.bam --minfasta 200000

# 七、分箱结果的物种及功能注释

## 7.1 分箱结果MAG的物种注释

### 7.1.1 使用GTDB-Tk进行MAG的基因组物种注释
    #MAG基因组分类及注释软件GTDB-Tk
    #使用conda创建新的环境并安装GTDB-Tk
    conda create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1 
    #激活GTDB-Tk软件安装环境
    conda activate gtdbtk-2.1.1 
    #直接使用软件自带脚本进行数据库配置，数据库配置位置为~/miniconda3/envs/gtdbtk-2.1.1/share/gtdbtk-2.1.1/db/gtdbtk_r207_v2_data.tar.gz
    download-db.sh

    #同时可选择手动配置数据库
    #设置数据库路径
    mkdir -p ~/db/gtdb & cd ~/db/gtdb
    #下载并解压数据库
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz  
    tar -zxvf gtdbtk_data.tar.gz
    #设置数据库位置，注意修改软件安装位置
    locate gtdbtk.sh # 查找配置文件位置
    #修改PATH=后面的路径为数据库解压目录，如/home/meta/db/gtdb/release95/
    vim /conda/envs/gtdbtk/etc/conda/activate.d/gtdbtk.sh

    #从官网下载适当的数据库，下载214版本的数据库
    wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz 
    #解压数据库
    tar -zcvf gtdbtk_r214_data.tar.gz
    #将下载的数据库配置到软件环境变量
    export GTDBTK_DATA_PATH=/ifs1/User/yongxin/db/gtdb 

    #确定软件需要的第三方程序以及数据库是否正常
    gtdbtk check_install 
    #运行bin物种分类及进化树构建
    gtdbtk classify_wf --genome_dir maxbin2_bins/ --extension fa  --skip_ani_screen --out_dir gtdbtk 
    gtdbtk convert_to_itol --input some_tree.tree --output itol.tree #转换进化树格式到itol

## 7.2 分箱结果MAG的功能注释
