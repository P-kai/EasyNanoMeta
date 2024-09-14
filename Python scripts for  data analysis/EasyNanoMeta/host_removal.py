import os
import argparse
import subprocess
import sys
# 创建host_removal_out文件夹并执行host_removal操作
def host_removal(host_removal_reference, num_threads):
    # 创建host_removal_out文件夹
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    if not os.path.exists(host_removal_out_folder):
        os.makedirs(host_removal_out_folder)
        print(f"Created folder: {host_removal_out_folder}")
    else:
        print(f"Folder already exists: {host_removal_out_folder}")

    # 获取参考文件的前缀名
    reference_prefix = os.path.splitext(os.path.basename(host_removal_reference))[0]

    # 步骤1：遍历 ./easynanometa_result/porechop_out 中的所有 fastq 文件并转换为 fasta 文件
    porechop_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "porechop_out")
    fasta_files = []  # 用于保存转换后的 fasta 文件路径

    for fastq_file in os.listdir(porechop_out_folder):
        if fastq_file.endswith("_output.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]  # 获取文件名前缀
            fastq_path = os.path.join(porechop_out_folder, fastq_file)
            # 生成输出的 fasta 文件路径，确保以 .fasta 结尾
            fasta_output_path = os.path.join(host_removal_out_folder, f"{file_prefix.replace('_output', '')}.fasta")
            fastq_to_fasta_cmd = (
                f"singularity exec easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit fq2fa {fastq_path} > {fasta_output_path}\""
            )
            try:
                subprocess.run(fastq_to_fasta_cmd, shell=True, check=True)
                print(f"Converted {fastq_path} to {fasta_output_path}")
                fasta_files.append(fasta_output_path)
            except subprocess.CalledProcessError as e:
               print(f"Error converting {fastq_path} to {fasta_output_path}: {e}")  
               
    # 步骤2：建立minimap2比对索引
    hostremoval_reference_directory = os.path.dirname(host_removal_reference)
    fasta_file_name = os.path.basename(host_removal_reference)
    minimap2_index_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -d {host_removal_out_folder}/{reference_prefix}.min /mnt/{fasta_file_name}\""
    try:
        print(f"Building minimap2 index for host genome {host_removal_reference}...")
        subprocess.run(minimap2_index_cmd, shell=True, check=True)
        print(f"Minimap2 index created: {host_removal_out_folder}/{reference_prefix}.min")
    except subprocess.CalledProcessError as e:
        print(f"Error during minimap2 index creation: {e}")
        return
    
    # 步骤3：遍历生成的 fasta 文件并使用minimap2进行比对
    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]

    # 使用minimap2进行比对
        minimap2_align_cmd = f"singularity exec easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -ax map-ont -t {num_threads} {host_removal_out_folder}/{reference_prefix}.min {fasta_file} -o {host_removal_out_folder}/{file_prefix}_minimap.sam\""
        try:
            print(f"Aligning {fasta_file} to host genome...")
            subprocess.run(minimap2_align_cmd, shell=True, check=True)
            print(f"Alignment completed: {host_removal_out_folder}/{file_prefix}_minimap.sam")
        except subprocess.CalledProcessError as e:
            print(f"Error during alignment: {e}")
            return

        # 步骤4：提取未比对到宿主的序列
        extract_unmaped_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools view -bS -@ {num_threads} -T /mnt/{fasta_file_name} -f 4 {host_removal_out_folder}/{file_prefix}_minimap.sam > {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam\""
        try:
            print(f"Extracting unmaped reads from {host_removal_out_folder}/{file_prefix}_minimap.sam...")
            subprocess.run(extract_unmaped_cmd, shell=True, check=True)
            print(f"Unmaped reads extracted: {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during unmaped reads extraction: {e}")
            continue

        # 步骤5：对bam文件进行排序
        sort_bam_cmd = f"singularity exec easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools sort -n {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam -o {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam\""
        try:
            print(f"Sorting BAM files...")
            subprocess.run(sort_bam_cmd, shell=True, check=True)
            print(f"Sorted BAM files: {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM file sorting: {e}")
            continue
        
        # 对bam文件进行转换为fastq
        bam_to_fastq_cmd = f"singularity exec easynanometa1_7.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && bedtools bamtofastq -i {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam -fq {host_removal_out_folder}/{file_prefix}_fitted_raw.fastq\""
        try:
            print(f"Converting BAM files to FASTQ files...")
            subprocess.run(bam_to_fastq_cmd, shell=True, check=True)
            print(f"Conversion completed: {host_removal_out_folder}/{file_prefix}_fitted_raw.fastq")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM to FASTQ conversion: {e}")
            continue
        
# 主函数
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
    parser.add_argument(
        "-host-removal-reference", "--host-removal-reference",
        required=True,
        help="Path to the reference host genome fasta file for host removal."
    )
    
    # 解析参数
    args = parser.parse_args()
    
    
    # 如果没有提供任何参数（即仅输入了脚本名），打印帮助信息
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)

    # 进行宿主去除
    host_removal(args.host_removal_reference, args.threads)
    
if __name__ == "__main__":
    main()
