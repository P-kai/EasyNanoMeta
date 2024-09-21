import os
import subprocess

# Flye functionality
def run_flye(input_dir, output_dir, singularity_image, flye_executable, threads=24):
    # Check if input directory exists
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    # Check if output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through all fastq files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # Get file name without extension
            output_path = os.path.join(output_dir, f"{sample_name}_flye")

            # Create output path if it does not exist
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # Generate Singularity Flye command
            command = [
                "singularity", "exec", "-B",
                f"{input_dir}:/mnt",
                singularity_image,
                flye_executable,
                "--nano-raw", f"/mnt/{file_name}",
                "--threads", str(threads),
                "--meta",
                "-o", output_path
            ]

            # Print the command to verify it
            print(f"Running Flye on {file_name} with output to {output_path}")
            print("Command:", " ".join(command))

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Flye assembly for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Flye on {file_name}: {e}")

# Kraken2 functionality
def run_kraken2(input_dir, output_dir, singularity_image, kraken2_executable, kraken2_db, threads=24):
    # Check if input directory exists
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    # Check if output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through all fastq files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # Get file name without extension

            # Generate output file paths
            report_file = os.path.join(output_dir, f"{sample_name}_kraken2_report")
            result_file = os.path.join(output_dir, f"{sample_name}_kraken2_result")

            # Generate Singularity Kraken2 command
            command = [
                "singularity", "exec", "-B",
                f"{kraken2_db}:/mnt", singularity_image,
                "/bin/bash", "-c",
                f"source /tools/Miniconda3/bin/activate base && "
                f"{kraken2_executable} --db /mnt/k2_standard/ "
                f"--report {report_file} "
                f"--output {result_file} "
                f"-t {threads} {input_file}"
            ]

            # Print the command to verify it
            print(f"Running Kraken2 on {file_name} with report output to {report_file}")
            print("Command:", " ".join(command))

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Kraken2 analysis for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Kraken2 on {file_name}: {e}")


# Abricate functionality
def run_abricate(input_dir, output_dir, singularity_image, abricate_executable, db='ncbi'):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fa'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]
            output_file = os.path.join(output_dir, f"{sample_name}_abricate_output.txt")

            command = [
                "singularity", "exec", "-B",
                f"{input_dir}:/mnt",
                singularity_image,
                "/bin/bash", "-c",
                f'"source /tools/Miniconda3/bin/activate abricate && '
                f'{abricate_executable} --db {db} /mnt/{file_name} > {output_file}"'
            ]
            print(f"Running Abricate on {file_name} with output to {output_file}")
            try:
                subprocess.run(" ".join(command), shell=True, check=True)
                print(f"Abricate analysis for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Abricate on {file_name}: {e}")
