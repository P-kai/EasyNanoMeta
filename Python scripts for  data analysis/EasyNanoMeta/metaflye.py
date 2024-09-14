import os
import subprocess
import argparse
import sys
import shutil

def metaflye(num_threads):
    # 创建metaflye_out文件夹
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    if not os.path.exists(metaflye_out_folder):
        os.makedirs(metaflye_out_folder)
        print(f"Created folder: {metaflye_out_folder}")
    else:
        print(f"Folder already exists: {metaflye_out_folder}")
    
    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
        
    # 遍历 host_removal_out 文件夹中的所有_fitted_raw.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            unique_fastq_path = os.path.join(host_removal_out_folder, f"{file_prefix}_unique.fastq")

            # 使用seqkit去除重复的读取ID
            seqkit_rmdup_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit rmdup {fastq_path} -o {unique_fastq_path}\""
            )
            try:
                print(f"Removing duplicates from {fastq_file}...")
                subprocess.run(seqkit_rmdup_cmd, shell=True, check=True)
                print(f"Duplicates removed. Output file: {unique_fastq_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error during seqkit processing for {fastq_file}: {e}")
                continue

    # 遍历 host_removal_out 文件夹中的所有_unique.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_unique.fastq"):
            file_prefix = fastq_file.replace("_unique.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            # 生成输出的文件路径
            metaflye_out_path = os.path.join(metaflye_out_folder, f"{file_prefix}_flye_out")

            # 步骤1：使用metaflye进行三代宏基因组批量组装
            metaflye_cmd = (
                f"singularity exec easynanometa1_10.sif /tools/Software/Flye-2.9.2/bin/flye "
                f"--nano-raw {fastq_path} "
                f"--threads {num_threads} "
                f"--meta -o {metaflye_out_path} "
            )
            try:
                print(f"Running metaFlye for {fastq_file}...")
                subprocess.run(metaflye_cmd, shell=True, check=True)
                print(f"metaFlye assembly completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during metaFlye assembly for {fastq_file}: {e}")
                continue

def mv_metaflye(metaflye_out_folder):
    # 遍历 metaflye_out_folder 中所有 _flye_out 文件夹
    for result_folder in os.listdir(metaflye_out_folder):
        if result_folder.endswith("_flye_out"):
            # 提取文件前缀，去除 '_flye_out'
            file_prefix = result_folder.replace("_flye_out", "")
            result_folder_path = os.path.join(metaflye_out_folder, result_folder)
            
            # 获取文件夹中的 assembly.fasta 文件路径
            assembly_fasta_path = os.path.join(result_folder_path, "assembly.fasta")
            
            if os.path.exists(assembly_fasta_path):
                # 生成新的 fasta 文件名，保存在与 assembly.fasta 相同的文件夹下
                new_fasta_name = f"{file_prefix}.fasta"
                new_fasta_path = os.path.join(result_folder_path, new_fasta_name)
                
                # 重命名 assembly.fasta 文件
                try:
                    shutil.move(assembly_fasta_path, new_fasta_path)
                    print(f"Renamed {assembly_fasta_path} to {new_fasta_path}")
                except Exception as e:
                    print(f"Error renaming {assembly_fasta_path}: {e}")
            else:
                print(f"assembly.fasta not found in {result_folder_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Find .fastq or .fq files in the specified folder and process them."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=24,
        help="Number of threads to use for all operations (default: 24)."
    )
    # 解析参数
    args = parser.parse_args()

    # 如果没有提供任何参数（即仅输入了脚本名），打印帮助信息
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)

    # 进行metaFlye分析
    metaflye(args.threads)
    
    # 移动并重命名 metaFlye 的输出文件
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    mv_metaflye(metaflye_out_folder)

if __name__ == "__main__":
    main()
