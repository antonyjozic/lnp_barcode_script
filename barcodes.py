#
# Author: Antony Jozic
#

from datetime import datetime
import gzip
import sys
import os

import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO


def get_fnames(root):
    """This function takes in a root directory path and walks through all subdirectories and files.
       The function checks each file in the subdirectories for the '.fastq.gz' file extension, adds then to a list and returns the list."""

    fnames = []
    for root, dirs, files in os.walk(root):
        for file in files:
            if file.endswith('.fastq.gz'):
                fnames.append(os.path.join(root, file))
    return fnames


def get_bc(barcode_path):
    """This function takes in the path to a barcode.xlsx file formatted in the same manner as the barcode.xlsx file found in this repository.
       The excel file is turned into a dataframe and then barcode sequences are extracted as a list.
       This list of barcode sequences is returned."""

    barcode_df = pd.read_excel(barcode_path)
    seq_bc = list(barcode_df['Barcode Seq:'])
    return seq_bc


def parse_fastq(fastq_path, bc_list):
    """This function takes in the path to a fastq.gz file and a list of barcodes.
       The function parses through each read and checks if the read contains any of the barcodes in the barcode list.
       If the read contains a barcode, the barcode count value for that barcode is iterated. A running tally of total reads is kept using a counter."""

    print("Starting analysis of: " + fastq_path + "\t" +
          datetime.today().strftime("%Y-%m-%d %H:%M:%S"))

    bc_counts = dict.fromkeys(bc_list, 0)
    with gzip.open(fastq_path, "rt") as f:
        for record in SeqIO.parse(f, "fastq"):
            for bc in bc_list:
                if bc in record:
                    bc_counts[bc] += 1
                    break
                else:
                    continue
    total_reads = sum(bc_counts.values())

    bc_counts["tot_read"] = total_reads
    fbc_df = pd.Series(bc_counts).to_frame(name="BC Read")

    print("Analysis of: " + fastq_path + " completed\t" +
          datetime.today().strftime("%Y-%m-%d %H:%M:%S"))
    return fbc_df


def main():
    """The main function pulls the paths of the 'fastq.gz' files given a root directory path. 
    The function also pulls the barcodes from an excel file and stores them in a list.
    Each 'fastq.gz' file is then analyzed and a final output file is created."""

    # pull parameters in
    fastq_fnames = get_fnames(sys.argv[1])
    bc_list = get_bc(sys.argv[2])
    rc = bool(int(sys.argv[3]))

    # name for output file
    df_fname = "barcode_reads.xlsx"

    # if reverse complement (revcomp) barcodes are required, revcomp the barcodes.xlsx file
    if rc == 1:
        rc_list = [str(Seq(x).reverse_complement()) for x in bc_list]
        bc_list = rc_list
        df_fname = "revcomp_barcode_reads.xlsx"

    # create final output file
    bc_df = pd.DataFrame()
    for fname in fastq_fnames:
        sample_bc_df = parse_fastq(fname, bc_list)
        colname = fname.split('/')[-1]
        bc_df[colname] = sample_bc_df
    bc_df.to_excel(df_fname)


usage = "To use this script, type: python3 <path to barcodes.py> <path to FASTQ.gz file directory> <path to barcode.xlsx file> <use reverse complement of barcodes (1 = yes / 0 = no)>"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()

main()
