# GP-transmission -- Influenza viral transmission bottlenecks and evolution
Welcome to the repository for the supporting information and files necessary to reproduce the analyses from Holmes et al (2025) "Viral expansion after transfer is a primary driver of influenza A virus transmission bottlenecks". 

The sequence data generated for this manuscript are available on the NCBI SRA under the accesssion number [PRJNA1253752](https://www.ncbi.nlm.nih.gov/bioproject/1253752).


## How this repository is organized

### barcode_analysis
This directory contains files that were generated during the sequence analysis of the barcoded virus used in this study, including the exact version of BarcodeID.py that was used for data processing (v1.2), and the script written to generate stacked bar plots shown in manuscript. 

The files are separated into analysis of samples collected from guinea pigs (barcode_info) and analysis of viral plaques isolated from the samples that were collected from guinea pigs (plaque_pick_barcode_info).

For a more detailed explanation of the output files generated from running BarcodeID.py, please visit the repository [here](https://github.com/Lowen-Lab/BarcodeID).

### WGS_analysis
This directory contains the files generated during analysis of whole genome sequence analysis of samples for our manuscript. The scripts used were a custom modified version of our lab's FluSAP pipeline (found [here](https://github.com/Lowen-Lab/FluSAP)). The files that are included here are those that are necessary to reproduce our within-host single nucleotide variant (iSNV) detection results, and to reproduce our analysis of WGS reads that span the barcoded region of the virus. 
