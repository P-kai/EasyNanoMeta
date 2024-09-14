import os
import argparse
import subprocess
import sys

# 定义 create_folder 函数
def create_folder(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
            print(f"Created folder: {folder}")
        except OSError as e:
            print(f"Error creating folder {folder}: {e}")
            sys.exit(1)
    else:
        print(f"Folder already exists: {folder}")
        
def abricate(num_threads):
 # 创建输出文件夹
    abricate_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "abricate_out")
    folders = {
        "arg": os.path.join(abricate_out_folder, "abricate_arg_out"),
        "vf": os.path.join(abricate_out_folder, "abricate_vf_out"),
        "is": os.path.join(abricate_out_folder, "abricate_is_out"),
        "abundance": os.path.join(abricate_out_folder, "arg_abundance_out")
    }
    for folder in folders.values():
        create_folder(folder)
            
    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")

    # 步骤1：遍历 ./easynanometa_result/host_removal_out 中的所有 fastq 文件并转换为 fasta 文件
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    fasta_files = []  # 用于保存转换后的 fasta 文件路径

    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]  # 获取文件名前缀
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            # 生成输出的 fasta 文件路径，确保以 .fasta 结尾
            fasta_output_path = os.path.join(abricate_out_folder, f"{file_prefix.replace('_fitted_raw', '')}.fasta")
            fastq_to_fasta_cmd = (
                f"singularity exec easynanometa1_9.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit fq2fa {fastq_path} > {fasta_output_path}\""
            )
            try:
                subprocess.run(fastq_to_fasta_cmd, shell=True, check=True)
                print(f"Converted {fastq_path} to {fasta_output_path}")
                fasta_files.append(fasta_output_path)
            except subprocess.CalledProcessError as e:
               print(f"Error converting {fastq_path} to {fasta_output_path}: {e}")  

    # 步骤2：获取每个 fasta 文件的基因组大小信息
    genome_sizes = {}  # 用于保存每个文件的基因组大小
    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
        # 使用 seqkit stats 获取基因组大小信息
        seqkit_cmd = (
            f"singularity exec easynanometa1_9.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
            f"seqkit stats {fasta_file} --all --tabular\""
        )
        try:
            result = subprocess.run(seqkit_cmd, shell=True, check=True, stdout=subprocess.PIPE, text=True)
            output = result.stdout.strip().split("\n")[-1]  # 获取最后一行的输出
            genome_size = output.split("\t")[1]  # 获取长度列的值
            genome_sizes[file_prefix] = genome_size
            print(f"Genome size of {fasta_file}: {genome_size}")
        except subprocess.CalledProcessError as e:
            print(f"Error getting genome size for {fasta_file}: {e}")
            return
            
    # 步骤2：遍历生成的 fasta 文件并使用abricate注释耐药基因
    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
                   
        # 使用abricate注释耐药基因
        abricate_arg_cmd = (
            f"singularity exec easynanometa1_9.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db ncbi -t {num_threads} {fasta_file} > {abricate_out_folder}/abricate_arg_out/{file_prefix}_ncbi_result\""          
        )          
        # 执行abricate_arg命令
        try:
            print(f"Annotatting ARGs of {fasta_file} ...")
            subprocess.run(abricate_arg_cmd, shell=True, check=True)
            print(f"Annotatting ARGs completed: {abricate_out_folder}/abricate_arg_out/{file_prefix}_ncbi_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during Annotatting ARGs: {e}")
            return
        
        # 步骤3：使用abricate注释毒力基因
        abricate_vf_cmd = (
            f"singularity exec easynanometa1_9.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db vfdb -t {num_threads} {fasta_file} > {abricate_out_folder}/abricate_vf_out/{file_prefix}_vf_result\""
        )
        # 执行abricate_vs命令
        try:
            print(f"Annotatting VF of {fasta_file} ...")
            subprocess.run(abricate_vf_cmd, shell=True, check=True)
            print(f"Annotatting VFs completed: {abricate_out_folder}/abricate_vf_out/{file_prefix}_vf_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during Annotatting VFs: {e}")
            return

        # 步骤4：使用abricate注释插入序列
        abricate_is_cmd = (
            f"singularity exec easynanometa1_9.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db ISfinder -t {num_threads} {fasta_file} > {abricate_out_folder}/abricate_is_out/{file_prefix}_is_result\""
        )
        # 执行abricate_is命令
        try:
            print(f"Annotatting IS of {fasta_file} ...")
            subprocess.run(abricate_is_cmd, shell=True, check=True)
            print(f"Annotatting IS completed: {abricate_out_folder}/abricate_is_out/{file_prefix}_is_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during Annotatting IS: {e}")
            return
        
        
    # 步骤5：遍历 abricate_arg_out 目录中的 {file_prefix}_ncbi_result 文件并计算丰度
    abricate_arg_out_folder = folders["arg"]
    arg_abundance_out_folder = folders["abundance"]
    for ncbi_result_file in os.listdir(abricate_arg_out_folder):
        if ncbi_result_file.endswith("_ncbi_result"):
            file_prefix = os.path.splitext(ncbi_result_file)[0]
            ncbi_result_path = os.path.join(abricate_arg_out_folder, ncbi_result_file)
            abundance_output_path = os.path.join(arg_abundance_out_folder, f"{file_prefix.replace('_ncbi_result', '')}_abundance_result")
            
            # 获取对应的基因组大小信息
            genome_size = genome_sizes.get(file_prefix, "unknown")
            if genome_size == "unknown":
                print(f"Genome size for {file_prefix} not found, skipping abundance calculation.")
                continue

            # 计算耐药基因丰度
            abundance_calculate_cmd = (
                f"singularity exec easynanometa1_9.sif /tools/Miniconda3/bin/python "
                f"/tools/Software/abundance_calculate.py --i {ncbi_result_path} --data_size {genome_size} > {abundance_output_path} "          
            )   
            try:
                print(f"Calculating abundance for {ncbi_result_file} with genome size {genome_size} ...")
                subprocess.run(abundance_calculate_cmd, shell=True, check=True)
                print(f"Abundance calculation completed: {abundance_output_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error during abundance calculation: {e}")
        
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
    
    # 解析参数
    args = parser.parse_args()
    
    
    # 如果没有提供任何参数（即仅输入了脚本名），打印帮助信息
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)
    
    # 进行abricate分析
    abricate(args.threads)

if __name__ == "__main__":
    main()
