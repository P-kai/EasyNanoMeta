import os
import argparse
import sys
import subprocess

def semi_bin(num_threads):
    # 创建输出文件夹
    semi_bin_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "semi_bin_out")
    bam_out_folder = os.path.join(semi_bin_out_folder, "bam_out")
    
    folders = {
        "bam": bam_out_folder
    }
    
    for folder in folders.values():
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
                print(f"Created folder: {folder}")
            except OSError as e:
                print(f"Error creating folder {folder}: {e}")
                sys.exit(1)
        else:
            print(f"Folder already exists: {folder}")

    # 固定 nextpolish_out_folder 路径为 ./easynanometa_result/nextpolish_out
    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")
    
    # 获取所有的前缀
    prefixes = []
    
    for folder_name in os.listdir(nextpolish_out_folder):
        folder_path = os.path.join(nextpolish_out_folder, folder_name)
        
        if os.path.isdir(folder_path):
            if folder_name.endswith("_nextpolish_out"):
                prefix = folder_name.rsplit("_nextpolish_out", 1)[0]
                prefixes.append((folder_path, prefix))
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
        
        for fasta_file in fasta_files:
            # 对组装基因组进行索引创建，获取排序的bam文件
            minimap2_index_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -d catalogue.mmi {fasta_file} \""
            )
            try:
                print(f"Creating index file for {fasta_file}...")
                subprocess.run(minimap2_index_cmd, shell=True, check=True)
                print(f"Index procedures completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during index procedures for {fasta_file}: {e}")
                continue

    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    # 遍历 host_removal_out 文件夹中的所有_fitted_raw.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.bam")
            
            # 比对获取bam文件
            minimap2_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -t {num_threads} -N 5 -ax map-ont catalogue.mmi {fastq_path} | samtools view -F 3584 -b --threads 8 > {bam_out_path} \""
            )
            try:
                print(f"BAM file generated for {fastq_path}...")
                subprocess.run(minimap2_cmd, shell=True, check=True)
                print(f"BAM file generated successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file generated procedures for {fastq_path}: {e}")
                continue
            
            # 对bam文件进行排序
            sorted_bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.sorted.bam")
            samtools_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"samtools sort -@ 10 {bam_out_path} > {sorted_bam_out_path} \""
            )
            try:
                print(f"BAM file sorted for {file_prefix}.bam...")
                subprocess.run(samtools_cmd, shell=True, check=True)
                print(f"BAM file sorted successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file sorted procedures for {file_prefix}.bam: {e}")
                continue
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
    
        for fasta_file in fasta_files:        
            # 使用SemiBin运行bin
            semi_bin_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate SemiBin && "
                f"SemiBin single_easy_bin -i {fasta_file}  --sequencing-type long_read -b {sorted_bam_out_path} -o {semi_bin_out_folder}/{file_prefix}_bin_out --environment global  \""
            )
            try:
                print(f"Perform bin for {fasta_file}...")
                subprocess.run(semi_bin_cmd, shell=True, check=True)
                print(f"Perform bin successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during bin procedures for {fasta_file}: {e}")
                continue

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
    args = parser.parse_args()
    
    semi_bin(args.threads)

if __name__ == "__main__":
    main()