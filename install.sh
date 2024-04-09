# 本脚本用于EasyNanoMeta分析流程的软件安装
# 设置数据库路径及软件安装路径

    db=~/db
    mkdir -p ${db} && cd ${db}
    # 软件安装位置Software installation location，默认为~/miniconda3，测试服务器为/anaconda3
    soft=~/miniconda3
    # 经常使用的服务器环境，可把全文${db}和${soft}替换为绝对路径，将不再需要每次读取以上环境变量
    # In the frequently used server environment, you can replace the variable ${db} and ${soft} with absolute paths, and you will no longer need to run the above environment variables every time
    # 可选：初始化环境变量，可能提高软件安装成功率
    # Optional: Initialize environment variables, which may improve the success rate of software installation
    PATH=${soft}/bin:${soft}/condabin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${db}/EasyMicrobiome/linux:${db}/EasyMicrobiome/script
    echo $PATH

# 下载最新版miniconda3，软件管理工具
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

# 安装数据质控及去宿主软件
## 安装porechop_abi
    # 创建porechop_abi单独的conda环境，使用conda安装软件
    conda create -y -n porechop_abi
    conda activate porechop_abi
    conda install -f -c conda-forge -c bioconda  porechop_abi

## 安装nanopack
    #创建nanopack环境,python3.10版本可成功安装，其他python版本可能安装失败。
    conda create -y -n nanopack python=3.10 
    conda activate nanopack
    
    #使用pip安装nanopack
    pip install nanopack

    #不成功，可尝试使用pip的国内清华源安装nanopack，由于网络问题的失败，可多尝试几次。
    pip install -i https://pypi.tuna.tsinghua.edu.cn/simple nanopack 

    #注意最终安装成功的软件路径提示，一般为path下
    /ifs1/User/yongxin/.local/bin/

## 安装minimap2，samtools，bedtools
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
