import os
import argparse
import subprocess
import sys
import shutil 
import pandas as pd

def gtdbtk(gtdbtk_db, num_threads):
    # 创建输出文件夹
    gtdbtk_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "gtdbtk_out")
    if not os.path.exists(gtdbtk_out_folder):
        os.makedirs(gtdbtk_out_folder)
        print(f"Created folder: {gtdbtk_out_folder}")
    else:
        print(f"Folder already exists: {gtdbtk_out_folder}")

    # 固定 semi_bin_out_folder 路径为 ./easynanometa_result/semi_bin_out
    semi_bin_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "semi_bin_out")
   
    # 获取所有的前缀
    prefixes = []
    
    for folder_name in os.listdir(semi_bin_out_folder):
        folder_path = os.path.join(semi_bin_out_folder, folder_name)
        
        if os.path.isdir(folder_path) and folder_name.endswith("_bin_out"):
            prefix = folder_name.rsplit("_bin_out", 1)[0]
            prefixes.append((folder_path, prefix))

            # 遍历 folder_path 中的所有 output_bins 文件夹
            for subfolder_name in os.listdir(folder_path):
                output_bins_folder = os.path.join(folder_path, subfolder_name)
                
                if os.path.isdir(output_bins_folder) and subfolder_name == "output_bins":
                    print(f"Processing {output_bins_folder} for prefix {prefix}")

                    # 执行 gtdbtk_cmd 
                    gtdbtk_cmd = (
                        f"singularity exec -B {gtdbtk_db}:/mnt easynanometa.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate gtdbtk-2.2.6 && export GTDBTK_DATA_PATH=/mnt && "
                        f"gtdbtk classify_wf --genome_dir {output_bins_folder} --extension fa --skip_ani_screen --out_dir {gtdbtk_out_folder}/{prefix}_gtdbtk_out\""
                    )
                    try:
                        print(f"Perform gtdbtk for {folder_path}...")
                        subprocess.run(gtdbtk_cmd, shell=True, check=True)
                        print(f"Perform gtdbtk successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during gtdbtk procedures for {folder_path}: {e}")
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
    parser.add_argument(
        "-gtdbtk-db", "--gtdbtk-database",
        type=str,
        help="Path to GTDB-Tk database path."
    )
    
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    gtdbtk(args.gtdbtk_database, args.threads)

if __name__ == "__main__":
    main()
