import os
import subprocess
import argparse

def checkm2(checkm2_db, num_threads):
    # 创建输出文件夹
    checkm2_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "checkm2_out")
    if not os.path.exists(checkm2_out_folder):
        os.makedirs(checkm2_out_folder)
        print(f"Created folder: {checkm2_out_folder}")
    else:
        print(f"Folder already exists: {checkm2_out_folder}")

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

                    checkm2_db_directory = os.path.dirname(checkm2_db)
                    db_file_name = os.path.basename(checkm2_db)

                    # 执行 checkm2_cmd 对 output_bins_folder
                    checkm2_cmd = (
                        f"singularity exec -B {checkm2_db_directory}:/mnt easynanometa1_10.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate checkm2 && "
                        f"checkm2 predict --database_path /mnt/{db_file_name} --threads {num_threads} --input {output_bins_folder}/* --output-directory {checkm2_out_folder}/{prefix}_checkm2_out\""
                    )
                    try:
                        print(f"Perform checkm2 for {folder_path}...")
                        subprocess.run(checkm2_cmd, shell=True, check=True)
                        print(f"Perform checkm2 successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during checkm2 procedures for {folder_path}: {e}")
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
        "-checkm2-db", "--checkm2-database",
        type=str,
        help="Path to checkm2 database path."
    )
    args = parser.parse_args()
    
    checkm2(args.checkm2_database, args.threads)

if __name__ == "__main__":
    main()
