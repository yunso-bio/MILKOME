#!/usr/bin/env python3

import os
import shutil
import logging

# Setup logging
logfile = 'process_bins.log'
logging.basicConfig(filename=logfile, level=logging.ERROR,
                    format='%(asctime)s %(levelname)s:%(message)s')
bins_path = "path to a list of bins"

for i in os.listdir("bins_for_das"):
    n = os.path.basename(i)
    summary_file = os.path.join("bins_for_das", i, f"{n}_DASTool_summary.tsv")
    
    try:
        with open(summary_file, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        logging.error(f"Error reading summary file {summary_file}: {e}")
        continue

    for line in lines:
        c1 = line.split()[0]
        c2 = line.split()[1]
        
        if c2 == "no-header_contig_bins.tsv":
            bin_name = f"{n}_SB_{c1}.fa.gz"
            src = os.path.join(bins_path, "semi_single_bin/bins", n, "output_bins", f"SemiBin_{c1}.fa.gz")
        else:
            bin_name = f"{n}_{c1}.fna"
            src = os.path.join(bins_path, "avamb/bins", n, "bins", f"{c1}.fna")
        
        dst = os.path.join("filtered_bins", bin_name)
        
        try:
            shutil.copy(src, dst)
        except Exception as e:
            logging.error(f"Error copying from {src} to {dst}: {e}")

print("finished")
