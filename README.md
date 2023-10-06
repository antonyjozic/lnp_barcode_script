# lnp_barcode_script
A script used to extract barcode sequences from barcoded LNP experiments in the Sahay Lab

## Running the script:
  - Inputs: 
    - path to the 'barcodes.py' file
    - path to the directory containing FASTQ files
    - path to 'barcodes.xlsx' file
    - whether to use the reverse complement of the barcode: "1 = yes, 0 = no"

## Output:
  - The script outputs an excel file with the barcode counts of each fastq file in the FASTQ directory listed by column, with each row linked to a unique barcode sequence found in the barcodes.xlsx
