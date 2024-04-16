# 本脚本用于EasyNanoMeta分析流程的软件安装
## 设置数据库路径及软件安装路径

    db=~/db
    mkdir -p ${db} && cd ${db}
    # 软件安装位置Software installation location，默认为~/miniconda3，测试服务器为/anaconda3
    soft=~/miniconda3
    # 其它软件安装位置Other softwares installation location，默认为~/tools
    tool=~/tools
    # 经常使用的服务器环境，可把全文${db}和${soft}替换为绝对路径，将不再需要每次读取以上环境变量
    # In the frequently used server environment, you can replace the variable ${db} and ${soft} with absolute paths, and you will no longer need to run the above environment variables every time
    # 可选：初始化环境变量，可能提高软件安装成功率
    # Optional: Initialize environment variables, which may improve the success rate of software installation
    PATH=${soft}/bin:${soft}/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${db}/EasyMicrobiome/linux:${db}/EasyMicrobiome/script
    echo $PATH

## 下载最新版miniconda3，软件管理工具
    wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    # 安装，-b批量，-f无提示，-p目录，许可协议打yes
    bash Miniconda3-latest-Linux-x86_64.sh -b -f -p ${soft}
    # 激活，然后关闭终端重开，提示符前出现(base)即成功
    ${soft}/condabin/conda init
    source ~/.bashrc
    # 查看版本，conda 23.7.3, python 3.11.4
    conda -V 
    python --version
    # 添加常用频道
    conda config --add channels bioconda # 生物软件
    conda config --add channels conda-forge # Highest priority

    # conda默认配置文件为 ~/.condarc 查看配置文件位置
    conda install mamba -c conda-forge -y
    mamba install pandas -c conda-forge -y
    mamba install conda-pack -c conda-forge -y
    conda config --set channel_priority strict
    conda config --show-sources
    # 查看虚拟环境列表 
    conda env list

# 1.安装数据质控及去宿主软件
## 1.1 长读数据质控软件
### 安装porechop_abi
    # 创建porechop_abi单独的conda环境，使用conda安装软件
    conda create -y -n porechop_abi
    conda activate porechop_abi
    conda install -f -c conda-forge -c bioconda  porechop_abi

### 安装nanopack
    #创建nanopack环境,python3.10版本可成功安装，其他python版本可能安装失败。
    conda create -y -n nanopack python=3.10 
    conda activate nanopack
    
    #使用pip安装nanopack
    pip install nanopack

    #不成功，可尝试使用pip的国内清华源安装nanopack，由于网络问题的失败，可多尝试几次。
    pip install -i https://pypi.tuna.tsinghua.edu.cn/simple nanopack 

    #注意最终安装成功的软件路径提示，一般为path下
    /ifs1/User/yongxin/.local/bin/

## 1.2 长读数据去宿主软件
### 安装minimap2，samtools，bedtools
    conda create -y -n host_removal
    conda activate host_removal
    conda install -c bioconda minimap2
    conda install -c bioconda samtools
    conda install -c bioconda bedtools

    #使用打包的conda package安装
    #package下载：https://figshare.com/account/projects/201156/articles/25569159
    mkdir ~/miniconda3/envs/host_removal/
    tar -xzvf host_removal.tar.gz -C ~/miniconda3/envs/host_removal/
    conda activate host_removal
    conda unpack

## 1.3 短读数据质控软件
### 安装fastp
    conda install fastp

## 1.4 数据格式转换：转换SRA数据到fastq格式
### 安装fastq-dump
    # 直接下载安装
    cd ~/tools
    #centos版本
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.6/sratoolkit.3.0.6-centos_linux64.tar.gz
    #ubuntu版本
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.6/sratoolkit.3.0.6-ubuntu64.tar.gz

    # 解压使用
    tar -zxvf sratoolkit.3.0.6-centos_linux64.tar.gz
    cd sratoolkit.3.0.6-centos_linux64

    # 查看软件版本：3.0.6
    ~/tools/sratoolkit.3.0.6-centos_linux64/bin/fastq-dump --version

# 2. 安装基于长读序列的物种注释、功能注释软件
## 2.1 长读序列物种注释
### 安装centrifuge
    #直接下载最新版本软件安装包，解压后编译
    cd ~/tools
    wget https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v1.0.4.tar.gz
    tar -zxvf v1.0.4.tar.gz
    cd centrifuge-1.0.4
    make
    make install prefix=~/tools/centrifuge-1.0.4

    #使用git克隆到本地进行编译安装
    git clone https://github.com/DaehwanKimLab/centrifuge
    cd centrifuge
    make
    make install prefix=~/tools/centrifuge-1.0.4

    #数据库配置，进入存储数据库的文件夹，进行数据库下载
    #h+p+v+c: human genome, prokaryotic genomes, and viral genomes including 106 SARS-CoV-2 complete genomes
    cd ~/db
    wget https://zenodo.org/record/3732127/files/h%2Bp%2Bv%2Bc.tar.gz?download=1
    tar -zxvf centrifuge_h+p+v.tar.gz

    #查看软件版本：1.0.4
    ~/tools/centrifuge-1.0.4/centrifuge --version

### 安装kraken2
    #直接下载软件安装包进行解压安装
    cd ~/tools
    wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.tar.gz
    tar -zxvf v2.1.3.tar.gz
    cd kraken2-2.1.3/
    sh install_kraken2.sh ~/tools/kraken2-2.1.3/

    #使用conda进行软件安装
    conda create -n kraken2
    conda activate kraken2
    conda install -y kraken2

    #数据库配置
    #直接使用kraken2自带脚本进行数据库下载
    kraken2-build --standard --threads 24 --db ~/db/kraken2_db

    #下载第三方构建的数据库直接使用，推荐网站：https://benlangmead.github.io/aws-indexes/k2
    #多种kraken2数据库可选，更新及时，并且可直接选择后下载解压使用，此处下载standard数据库
    cd ~/db
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz

    #解压数据库到指定数据库位置
    tar -zcvf k2_standard_20230605.tar.gz -C ~/db/k2_standard/

    #查看软件版本:2.1.3
    ~/tools/kraken2-2.1.3/bin/kraken2 --version

## 2.2 长读序列功能注释
### 安装abricate
#直接使用conda安装abricate
    #创建abricate单独的环境
    conda create -n abricate
    #激活abricate环境并进行软件安装，指定安装1.0.1版本，默认安装会安装低版本
    conda activate abricate
    conda install -y -c bioconda abricate=1.0.1

    #查看abricate当前的数据库
    abricate --list

    #软件提供card和ncbi AMRFinder两个抗性数据库，首先更新已有的数据库
    abricate-get_db --db ncbi #ncbi表示AMRFinder数据库
    abricate-get_db --db card

    #此外，也可自己下载相应的数据库，进行构建
    #进入abricate数据库目录，若为conda或miniconda安装，则在conda目录下
    cd /path/to/abricate/db
    #创建数据库文件夹，并将fasta格式的数据库拷贝到文件夹内，重命名为sequences
    mkdir your_database_name
    cp /your/database/database.fasta your_database_name/sequences
    #创建数据库，名称为 your_database_name
    makeblastdb -in sequences -title your_database_name -dbtype nucl -hash_index

# 3. 安装长读、混合宏基因组组装软件
## 3.1 安装长读宏基因组组装软件
### 安装MetaFlye
    # 软件下载
    cd ~/tools
    wget https://github.com/fenderglass/Flye/archive/refs/tags/2.9.2.tar.gz
    # 解压使用
    tar -zxvf Flye-2.9.2.tar.gz
    # 查看软件版本，版本：2.9.2-b1786
    ~/tools/Flye-2.9.2/bin/flye --version

### 安装Canu
    # 软件下载，链接：https://github.com/marbl/canu/releases
    cd ~/tools
    curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz
    #解压软件，即可直接使用
    #进入软件安装及执行目录，并查看软件安装目录和版本，版本：canu 2.2
    tar -xJf canu-2.2.Linux-amd64.tar.xz 
    cd canu-2.2/bin && pwd 
    ~/tools/canu-2.2/bin/canu
    ~/tools/canu-2.2/bin/canu --version

### 安装wdtbg2
    # 软件下载及安装
    cd ~/tools
    git clone https://github.com/ruanjue/wtdbg2
    cd wtdbg2 && make

    # 安装wtdbg2依赖工具：samtools，minimap2
    conda install -y samtools
    conda install -y minimap2

    # 查看软件版本：wtdbg2 2.5
    ~/tools/wtdbg2/wtdbg2 --version

### 安装NextDenovo
    #软件下载及安装
    wget https://github.com/Nextomics/NextDenovo/releases/latest/download/NextDenovo.tgz 
    tar -vxzf NextDenovo.tgz

    #安装软件依赖包
    pip install paralleltask

    #查看软件版本：nextDenovo 2.5.2
    ~/tools/NextDenovo/nextDenovo --version
    
## 3.2 安装长短读混合宏基因组组装软件
### 安装OPERA-MS
    #软件安装
    #使用conda配置软件安装单独环境，安装软件依赖的perl模块
    conda create -n operams python=3.9
    #激活conda环境 
    conda activate operams
    #在conda环境中安装依赖的perl模块
    conda install -c conda-forge perl-app-cpanminus 
    conda install -c compbiocore perl-switch perl==5.26.2
    conda install -c bioconda perl-file-which perl-statistics-basic perl-statistics-r

    #下载软件安装包
    cd ~/tools
    git clone https://github.com/CSB5/OPERA-MS.git
    #进入软件目录并执行软件编译，最后检查软件依赖的所有perl模块
    cd OPERA-MS
    make
    perl OPERA-MS.pl check-dependency

    #可能遇到问题 “Can't locate Switch.pm”
    #解决：寻找当前用户目录下有没有Switch.pm模块的安装 find ~/ -name "Switch.pm" 
    #将找到的模块写入perl路径中，例如： export PERL5LIB=~/perl5/lib/perl5/

    #配置OPERA-MS软件数据库
    perl OPERA-MS.pl install-db 

    #激活operams环境，查看软件版本：OPERA-MS v0.9.0
    conda activate operams
    perl ~/tools/OPERA-MS/OPERA-MS.pl

### 安装MetaSPAdes
    # 下载预编译的软件安装包，解压
    cd ~/tools
    wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz
    tar -zxvf SPAdes-3.15.5-Linux.tar.gz

    #查看软件版本：SPAdes genome assembler v3.15.5
    ~/tools/SPAdes-3.15.5-Linux/bin/spades.py --version

### 安装MetaPlatanus
    #使用conda进行软件安装
    #使用conda创建软件安装单独环境，进入环境进行软件安装
    conda create -n metaplatanus -y
    conda activate metaplatanus
    conda install -c conda-forge -c bioconda metaplatanus

    #查看软件版本：metaplatanus version v1.3.1
    metaplatanus --version

    #使用conda的package进行软件安装
    #package下载：https://figshare.com/account/projects/201156/articles/25573755
    cd ~/tools
    wget -c --no-check-certificate --no-proxy https://figshare.com/ndownloader/files/45563190 -O metaplatanus.tar.gz
    mkdir ~/miniconda3/envs/metaplatanus/
    tar -xzvf metaplatanus.tar.gz -C ~/miniconda3/envs/metaplatanus/
    conda activate metaplatanus
    conda unpack

### 安装Unicycler
    #使用conda进行软件安装
    #使用conda创建软件安装单独环境，进入环境进行软件安装
    onda create -n unicycler -y
    conda activate unicycler
    conda install unicycler -c bioconda -y

    #查看软件版本：Unicycler v0.5.0
    unicycler --version

    #使用conda的package进行软件安装
    #package下载：https://figshare.com/articles/software/Untitled_Item/25573875
    cd ~/tools
    wget -c --no-check-certificate --no-proxy https://figshare.com/ndownloader/files/45563385 -O unicycler.tar.gz
    mkdir ~/miniconda3/envs/unicycler/
    tar -xzvf unicycler.tar.gz -C ~/miniconda3/envs/unicycler/
    conda activate unicycler
    conda unpack

# 4. 安装组装结果校准软件
## 4.1 三代组装结果校准
### 安装minimap2和recon
    #使用minimap2与racon进行组装基因组校准，若当前环境无minimap2与recon，则首先使用conda进行软件安装
    #使用conda安装minimap2与recon，已安装软件可跳过此步骤
    conda install -y minimap2
    conda install -y racon

    #查看软件版本：minimap2 2.17-r943-dirty；
    minimap2 --version

### 安装NextPolish
    #下载软件安装包，解压并编译
    cd ~/tools
    wget https://github.com/Nextomics/NextPolish/releases/latest/download/NextPolish.tgz
    tar -vxzf NextPolish.tgz && cd NextPolish && make

    #安装软件依赖
    pip install paralleltask

    #查看软件版本：v1.4.1
    ~/tools/NextPolish/nextPolish --version

## 4.2 二代或者二三代混合组装结果校准
### 安装Pilon
    # 直接下载，然后使用java调用
    cd ~/tools
    wget https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar

    # 查看pilon版本：Pilon version 1.24 Thu Jan 28 13:00:45 2021 -0500
    java -Xmx16G -jar ~/tools/pilon-1.24.jar --version

### 安装bwa
    # 使用git克隆bwa后编译
    git clone https://github.com/lh3/bwa.git
    cd bwa; make
    
    #添加环境变量
    export PATH=~/tools/bwa:$PATH

# 5. 安装宏基因组分箱软件
## 5.1 纳米孔宏基因组分箱软件安装
### 安装SemiBin
    #使用conda创建单独的环境进行软件安装
    conda create -n SemiBin
    conda activate SemiBin
    conda install -c conda-forge -c bioconda semibin

    #查看软件版本：2.1.0
    SemiBin2 --version

    #使用conda的package进行软件安装
    #package下载：https://figshare.com/account/projects/201156/articles/25574010
    cd ~/tools
    wget -c --no-check-certificate --no-proxy https://figshare.com/ndownloader/files/45563634 -O SemiBin.tar.gz
    mkdir ~/miniconda3/envs/SemiBin/
    tar -xzvf SemiBin.tar.gz -C ~/miniconda3/envs/SemiBin/
    conda activate SemiBin
    conda unpack

### 安装vamb
    #使用pip安装vamb
    pip install vamb -i https://pypi.tuna.tsinghua.edu.cn/simple
    
    #查看软件版本：Vamb 4.1.3
    vamb --version

### 安装metawrap
    #使用conda创建单独的环境进行软件安装
    conda create -y -n metawrap-env python=2.7
    conda activate metawrap-env
    
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels ursky
    
    conda install -y mamba 
    mamba install --only-deps -c ursky metawrap-mg

    #查看软件版本：1.3.2
    metawrap --version

    #使用conda的package进行软件安装
    #package下载：https://figshare.com/articles/software/MetaWrap_1_3_2/25603155
    cd ~/tools
    wget -c --no-check-certificate --no-proxy https://figshare.com/ndownloader/files/45651492 -O metawrap.tar.gz
    mkdir ~/miniconda3/envs/metawrap/
    tar -xzvf metawrap.tar.gz -C ~/miniconda3/envs/metawrap/
    conda activate metawrap
    conda unpack

# 6. 安装MAGs质控、物种注释、功能注释软件
## 6.1 MAGs质控软件
### 安装checkm
    # checkm已被整合到流程metawrap中，可直接激活metawrap环境进行使用
    # 也可使用conda创建环境单独安装checkm
    # 创建python=3.9的conda环境
    conda create -n checkm python=3.9
    conda activate checkm

    # 使用pip3安装checkm及其依赖环境
    pip3 install numpy
    pip3 install matplotlib
    pip3 install pysam
    pip3 install checkm-genome

    # 查看软件版本：v1.2.2
    checkm

### 安装checkm2
    # 使用conda安装checkm2
    # 创建checkm2环境，注意python版本为3.8，否则可能安装失败
    conda create -n checkm2 python=3.8
    conda activate checkm2
    # 使用mamba安装
    mamba install -c bioconda -c conda-forge checkm2

    # 使用conda的package进行软件安装
    #package下载：
    cd ~/tools
    wget -c --no-check-certificate --no-proxy  -O checkm2.tar.gz
    mkdir ~/miniconda3/envs/checkm2/
    tar -xzvf checkm2.tar.gz -C ~/miniconda3/envs/checkm2/
    conda activate checkm2
    conda unpack

    # 配置checkm2数据库
    # 直接使用checkm2脚本进行数据库下载（经常失败，推荐使用wget）
    checkm2 database --download
    
    # 使用wget进行数据库下载
    mkdir ~/db/checkm2 && cd ~/db/checkm2
    wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz
    # 解压数据库
    tar -zxvf checkm2_database.tar.gz

## 6.2 MAGs物种注释
### 安装gtdbtk
    # 使用conda进行软件安装
    conda create -n gtdbtk-2.2.6 -c conda-forge -c bioconda gtdbtk=2.2.6
    # 查看软件版本: v2.2.6
    conda activate gtdbtk-2.2.6
    gtdbtk --version

    # 使用conda的package进行软件安装
    #package下载：
    cd ~/tools
    wget -c --no-check-certificate --no-proxy  -O checkm2.tar.gz
    mkdir ~/miniconda3/envs/gtdbtk-2.2.6/
    tar -xzvf gtdbtk-2.2.6.tar.gz -C ~/miniconda3/envs/gtdbtk-2.2.6/
    conda activate gtdbtk-2.2.6
    conda unpack

    # 配置软件数据库
    mkdir -p ~/db/gtdbtk && cd ~/db/gtdbtk
    wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz
    tar -zxvf gtdbtk_data.tar.gz

## 6.3 MAGs功能注释
### 安装prokka
    # 使用conda安装prokka
    conda create -n prokka
    conda install -c conda-forge -c bioconda -c defaults prokka
