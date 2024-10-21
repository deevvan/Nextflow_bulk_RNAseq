# python srp_downloader.py --output_dir /path/to/output 
# --metadata_file /path/to/metadata.csv

import os
import multiprocessing
import pandas as pd
import argparse

def download_sra_run(srr_id, output_dir):
    # Download the file using fasterq-dump
    os.system(f"fasterq-dump -O {output_dir} {srr_id}")
    print(f"{srr_id} downloaded.")

def download_srr_files(output_dir, metadata_file):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the metadata file
    metadata_df = pd.read_csv(metadata_file)

    # Get the list of SRR IDs from metadata_df
    srr_ids = metadata_df['Run'].tolist()

    # Set the number of processes to use (adjust according to your CPU cores)
    num_processes = multiprocessing.cpu_count()

    # Create a multiprocessing pool
    pool = multiprocessing.Pool(processes=num_processes)

    # Map the download function to the SRR IDs using the multiprocessing pool
    pool.starmap(download_sra_run, [(srr_id, output_dir) for srr_id in srr_ids])

    # Close the pool
    pool.close()
    pool.join()

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Download SRR files using fasterq-dump')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the downloaded SRR files')
    parser.add_argument('-p', '--metadata_file', required=True, help='CSV file containing metadata file with SRR IDs')
    
    args = parser.parse_args()

    # Call the download function with the parsed arguments
    download_srr_files(args.output_dir, args.metadata_file)