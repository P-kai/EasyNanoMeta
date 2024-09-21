#!/usr/bin/env python

# easynanometa2.py

import argparse
from dep import run_flye, run_kraken2, run_abricate  # Import functions from dep.py

# Main function to choose Flye or Kraken2 or abricate
def main():
    # Create a top-level parser showing flye, kraken2 and abricate options
    parser = argparse.ArgumentParser(prog='easynano', description='Batch execute Flye, Kraken2 or abricate for multiple fastq files.')
    subparsers = parser.add_subparsers(dest='tool', help='Choose the tool to run.')

    # Create a subcommand parser for flye
    parser_flye = subparsers.add_parser('flye', help='Run Flye assembler')
    parser_flye.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_flye.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Flye output.')
    parser_flye.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_flye.add_argument('-f', '--flye_executable', default='/tools/Software/Flye-2.9.2/bin/flye', help='Path to Flye executable within the Singularity image.')
    parser_flye.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # Create a subcommand parser for kraken2
    parser_kraken2 = subparsers.add_parser('kraken2', help='Run Kraken2 classifier')
    parser_kraken2.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_kraken2.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Kraken2 output.')
    parser_kraken2.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_kraken2.add_argument('-k', '--kraken2_executable', default='/tools/Software/kraken2-2.1.3/kraken2', help='Path to Kraken2 executable within the Singularity image.')
    parser_kraken2.add_argument('-db', '--kraken2_db', default='/backup/database/kraken2', help='Path to the Kraken2 database.')
    parser_kraken2.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # Create a subcommand parser for abricate
    parser_abricate = subparsers.add_parser('abricate', help='Run abricate to identify functional genes')
    parser_abricate.add_argument('-i', '--input_dir', required=True, help='Directory containing the input .fa files.')
    parser_abricate.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Abricate output.')
    parser_abricate.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_abricate.add_argument('-e', '--abricate_executable', default='abricate', help='Path to Abricate executable within the Singularity image.')
    parser_abricate.add_argument('-d', '--db', default='ncbi', help='The database to use for Abricate (default: ncbi).')

    args = parser.parse_args()

    # Choose tool
    if args.tool == 'flye':
        run_flye(args.input_dir, args.output_dir, args.singularity_image, args.flye_executable, args.threads)
    elif args.tool == 'kraken2':
        run_kraken2(args.input_dir, args.output_dir, args.singularity_image, args.kraken2_executable, args.kraken2_db, args.threads)
    elif args.tool == 'abricate':
        run_abricate(args.input_dir, args.output_dir, args.singularity_image, args.abricate_executable, args.db)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
