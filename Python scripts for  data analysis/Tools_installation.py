import os
import argparse
import subprocess

def run_command(command, error_message):
    result = subprocess.run(command, capture_output=True, text=True, shell=True)
    if result.returncode != 0:
        raise Exception(f"{error_message}: {result.stderr.strip()}")
    return result.stdout.strip()

def Installation(checkm2_database, gtdbtk_database):
    try:
        # 检查是否已经安装了所需的软件包
        print("Checking installations...")

        # 检查git是否安装
        run_command('git --version', "Git is not installed. Please install Git first.")
        print("Git is installed.")

        # 检查Python版本
        python_version = run_command('python --version', "Python is not installed. Please install Python first.")
        print(f"Python version: {python_version}")

        # 检查conda是否安装
        run_command('conda --version', "Conda is not installed. Please install Conda first.")
        print("Conda is installed.")

        # 创建conda环境
        print("Creating conda environment 'pipeline1'...")
        run_command('conda create -n pipeline1 -y', "Failed to create conda environment.")
        print("Conda environment 'pipeline1' created successfully.")

        # 安装Flye软件包
        print("Installing Flye version 2.9.2 in 'pipeline1' conda environment...")
        run_command('conda install flye=2.9.2 -n pipeline1 -y', "Failed to install Flye.")
        print("Flye version 2.9.2 installed successfully in 'pipeline1' conda environment.")

        # 安装NextPolish软件包
        print("Installing NextPolish in 'pipeline1' conda environment...")
        run_command('conda install bioconda::nextpolish -n pipeline1 -y', "Failed to install NextPolish.")
        print("NextPolish installed successfully in 'pipeline1' conda environment.")

        # 安装Minimap2软件包
        print("Installing Minimap2 in 'pipeline1' conda environment...")
        run_command('conda install bioconda::minimap2 -n pipeline1 -y', "Failed to install Minimap2.")
        print("Minimap2 installed successfully in 'pipeline1' conda environment.")

        # 安装SemiBin软件包
        print("Installing SemiBin version 2.1.0 in 'pipeline1' conda environment...")
        run_command('conda install -c conda-forge -c bioconda semibin=2.1.0 -n pipeline1 -y', "Failed to install SemiBin.")
        print("SemiBin version 2.1.0 installed successfully in 'pipeline1' conda environment.")

        # 安装gtdbtk软件包
        print("Installing gtdbtk version 2.2.6 in 'pipeline1' conda environment...")
        run_command('conda install -c conda-forge -c bioconda gtdbtk=2.2.6 -n pipeline1 -y', "Failed to install gtdbtk.")
        print("gtdbtk version 2.2.6 installed successfully in 'pipeline1' conda environment.")

        # 配置gtdbtk数据库
        print(f"Configuring gtdbtk database at {gtdbtk_database}...")
        run_command(f'mkdir -p {gtdbtk_database} && cd {gtdbtk_database}', "Failed to create gtdbtk database directory.")
        run_command(f'wget -c https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_data.tar.gz -O {gtdbtk_database}/gtdbtk_data.tar.gz', "Failed to download gtdbtk database.")
        run_command(f'tar -zxvf {gtdbtk_database}/gtdbtk_data.tar.gz -C {gtdbtk_database}', "Failed to decompress gtdbtk database.")
        print("gtdbtk database downloaded and configured successfully.")

        # 安装CheckM2软件包
        print("Installing CheckM2 in 'pipeline1' conda environment...")
        run_command('mamba install -c bioconda -c conda-forge checkm2 -n pipeline1 -y', "Failed to install CheckM2.")
        print("CheckM2 installed successfully in 'pipeline1' conda environment.")

        # 下载并配置CheckM2数据库
        print(f"Downloading and configuring CheckM2 database at {checkm2_database}...")
        run_command(f'mkdir -p {checkm2_database} && cd {checkm2_database}', "Failed to create database directory.")
        run_command(f'wget https://zenodo.org/record/5571251/files/checkm2_database.tar.gz -O {checkm2_database}/checkm2_database.tar.gz', "Failed to download CheckM2 database.")
        run_command(f'tar -zxvf {checkm2_database}/checkm2_database.tar.gz -C {checkm2_database}', "Failed to decompress CheckM2 database.")
        print("CheckM2 database downloaded and configured successfully.")

        print("All checks and installations completed.")
        
    except Exception as e:
        print(f"An error occurred during installation checks: {e}")

def main():
    parser = argparse.ArgumentParser(description="Installation Check Script")
    parser.add_argument('--check', action='store_true', help="Run installation checks and installations")
    parser.add_argument('--checkm2_database', type=str, required=True, help="Absolute path to CheckM2 database directory")
    parser.add_argument('--gtdbtk_database', type=str, required=True, help="Absolute path to gtdbtk database directory")
    args = parser.parse_args()

    if args.check:
        Installation(args.checkm2_database, args.gtdbtk_database)

if __name__ == "__main__":
    main()
