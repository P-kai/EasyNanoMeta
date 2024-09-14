import os
import subprocess
import argparse
import shutil

def nextpolish(num_threads):
    # 创建nextpolish_out文件夹
    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")
    if not os.path.exists(nextpolish_out_folder):
        os.makedirs(nextpolish_out_folder)
        print(f"Created folder: {nextpolish_out_folder}")
    else:
        print(f"Folder already exists: {nextpolish_out_folder}")
    
    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    # 固定 metaflye_out_folder 路径为 ./easynanometa_result/metaflye_out
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    
    # 遍历 host_removal_out 文件夹中的所有_fitted_raw.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            lgs_out_path = os.path.join(nextpolish_out_folder, f"{file_prefix}.fofn")

            # 生成.fofn文件
            try:
                with open(lgs_out_path, 'w') as lgs_file:
                    lgs_file.write(fastq_path + "\n")
                print(f"fofn file created successfully for {fastq_file}. Output file: {lgs_out_path}")
            except Exception as e:
                print(f"Error during fofn file creation for {fastq_file}: {e}")
                continue

            # 生成.cfg文件
            cfg_file_path = os.path.join(nextpolish_out_folder, f"{file_prefix}.cfg")
            metaflye_prefix_path = os.path.join(metaflye_out_folder, f"{file_prefix}_flye_out")
            genome_path = os.path.join(metaflye_prefix_path, f"{file_prefix}.fasta")
            workdir_path = os.path.join(nextpolish_out_folder, f"{file_prefix}_nextpolish_out")
            cfg_content = f"""\

[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = {genome_path} #组装结果文件
genome_size = auto
workdir = {workdir_path}
polish_options = -p {num_threads}

[lgs_option]
lgs_fofn = {lgs_out_path}
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont
"""

            try:
                with open(cfg_file_path, 'w') as cfg_file:
                    cfg_file.write(cfg_content)
                print(f"cfg file created successfully. Output file: {cfg_file_path}")
            except Exception as e:
                print(f"Error during cfg file creation for {lgs_out_path}: {e}")

            # 执行nextpolish校准            
            nextpolish_cmd = (
                f"singularity exec easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"/tools/Software/NextPolish/nextPolish {cfg_file_path}\""
            )
            
            try:
                print(f"Performing calibration procedures for {file_prefix}.fasta...")
                subprocess.run(nextpolish_cmd, shell=True, check=True)
                print(f"Calibration procedures completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during calibration procedures for {file_prefix}.fasta: {e}")
                continue

def mv_nextpolish(nextpolish_out_folder):
    # 遍历 nextpolish_out_folder 中所有 _nextpolish_out 文件夹
    for result_folder in os.listdir(nextpolish_out_folder):
        if result_folder.endswith("_nextpolish_out"):
            # 提取文件前缀，去除 '_nextpolish_out'
            file_prefix = result_folder.replace("_nextpolish_out", "")
            result_folder_path = os.path.join(nextpolish_out_folder, result_folder)
            
            # 获取文件夹中的 assembly.fasta 文件路径
            calibration_fasta_path = os.path.join(result_folder_path, "genome.nextpolish.fasta")
            
            if os.path.exists(calibration_fasta_path):
                # 生成新的 fasta 文件名，保存在与 genome.nextpolish.fasta 相同的文件夹下
                new_fasta_name = f"{file_prefix}_nextpolish.fasta"
                new_fasta_path = os.path.join(result_folder_path, new_fasta_name)
                
                # 重命名 genome.nextpolish.fasta 文件
                try:
                    shutil.move(calibration_fasta_path, new_fasta_path)
                    print(f"Renamed {calibration_fasta_path} to {new_fasta_path}")
                except Exception as e:
                    print(f"Error renaming {calibration_fasta_path}: {e}")
            else:
                print(f"genome.nextpolish.fasta not found in {result_folder_path}")

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
    
    nextpolish(args.threads)

    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")
    
    mv_nextpolish(nextpolish_out_folder)

if __name__ == "__main__":
    main()
