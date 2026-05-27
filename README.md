# ORI-marker-screen

Analyses and data for identifying and quantifying abundances of origins of replication in a pooled ORI-marker screen

## Overview

This repository contains code and data for reproducing the analyses described in "A scalable framework for high-throughput identification of functional origins of replication in non-model bacteria". It also contains scripts and examples for how to identify and quantify the origins of replications present in a pool of plasmids, either after or before a selective screen. 

Conducting a plasmid pool screen and having issues with these scripts? Open an issue and we will respond to help ASAP!

## Updated branch: for ingest. 

Example commands:
```
python barcode_quantification.py -m test_data/test2.csv -d ./test_data/test2/
```

```
python barcode_quantification.py -m test_data/test1.csv -d ./test_data/amplicon_fastq
```

Output:

```
barcode_results.tsv - raw counts for the run
barcode_stats.tsv - statistics on the run
portal_ingest.tsv - table ready for ingest 
```

## Requirements

1. Python library requirements are listed in `requirements.txt`.
2. [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (both amplicon and whole plasmid sequencing)
3. [Bowtie2](https://github.com/BenLangmead/bowtie2) (whole plasmid sequencing)
4. [BBmerge](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) (amplicon sequencing)
5. [VSEARCH](https://github.com/torognes/vsearch) (amplicon sequencing)

## Install and quick start

```
git clone https://github.com/cultivarium/ORI-marker-screen.git
pip install -r requirements.txt
python barcode_quantification.py -d ms_data/amplicon_fastq/ -m ./ms_data/barcode_mapping_test.csv
```

## Data

Data files specific for the Cultivarium Possum Toolkit library of ORIs are available as follows:

- `barcode_references.fasta` - Barcode sequences associated with each ORI in the library in FASTA format.
- `barcodes.csv` - A list of the Barcodes associated with each ORI in the library.
- `origins.tsv` - Start and stop locations of each ORI region on each plasmid in the library.
- `./pool_data/*` - FASTA and Genbank files describing the plasmids (and their respective ORIs) within the Cultivarium Possum Toolkit. 
- Raw sequencing data for our preprint can be obtained through AWS S3: `aws s3 cp --recursive s3://cultivarium-sequencing/ORI-MARKER-RAW-DATA-MAY2023/ .`

## Amplicon sequencing for barcoded ORI quantification

For identifying the presence of amplicon barcodes from the Cultivarium Possum Toolkit. Each barcode is linked to an ORI within the library, and therefore the barcodes are used to identify the presence and abundances of each ORI within the pool. 

Example usage: `python barcode_quantification.py -d ms_data/amplicon_fastq/ -m ./ms_data/barcode_mapping_test.csv`.

All arguments:

```
usage: barcode_quantification.py [-h] -d FASTQ_DIRECTORY -m MAPPING_FILE [-l LIBRARY_INFO] [-b BBMAP_FOLDER] [-o OUTPUT_FOLDER]
                                 [--unmerged_reads] [--nanopore]

Amplicon BarSeq plasmid sequencing of a plasmid ORI pool.

options:
  -h, --help            show this help message and exit
  -d FASTQ_DIRECTORY, --fastq_directory FASTQ_DIRECTORY
                        Directory of FASTQ files. File names must take the form: sample_*_R1_*.fastq.gz
                        (or sample*.fastq[.gz], one file per sample, with --nanopore)
  -m MAPPING_FILE, --mapping_file MAPPING_FILE
                        Mapping file of comma separated columns FileName,Sample,Strain,Library.
  -l LIBRARY_INFO, --library_info LIBRARY_INFO
                        Library info file of comma separated information about ORIs in the library. If running input samples, 'Negative control
                        cutoff' can be blank or ignored. 
  -b BBMAP_FOLDER, --bbmap_folder BBMAP_FOLDER
                        Directory containing BBTools on your system
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Directory to store output files
  --unmerged_reads      Also process unmerged R1 (useful with 2x75 bp or lower quality reads)
  --nanopore            Single-end Nanopore amplicon mode (see below)
```

### Nanopore amplicon reads (`--nanopore`)

The same barcode amplicons can be sequenced single-end on Oxford Nanopore. Because
each Nanopore read already spans the entire amplicon, no read pairing/merging is
needed, and the reads are not strand-oriented. Pass `--nanopore` to handle this:

- **Input:** one FASTQ per sample (gzipped or not), matched as `FileName*.fastq[.gz]`
  (instead of paired `*_R1_*.fastq.gz` / `*_R2_*.fastq.gz`).
- **Trimming:** BBduk runs single-end with the same adapter/quality/entropy
  parameters (the R2-only `forcetrimright2` is dropped).
- **Merging:** skipped — the cleaned reads are taken directly as the full-length
  amplicons (they are converted to FASTA with `reformat.sh`).
- **Matching:** VSEARCH searches **both strands** (`--strand both`), since Nanopore
  reads arrive in both orientations. The identity threshold is unchanged (0.95).

All outputs are identical in format to the Illumina mode. Example:

```
python barcode_quantification.py -d ./nanopore_fastq -m mapping.csv -l library_info.csv --nanopore
```

Example outputs:

`./ms_data/barcode_stats.tsv` - Provides statistics on read filtering, merging, and matching for each sample. Have a look at this file to understand general quality of your run and identify any potential issues.

`./ms_data/barcode_results.tsv` - The read pair / amplicon counts for each ORI within each sample.

## Running on input libraries

The pipeline can now be run on input plasmid libraries, not just strains. To do so:

1. first, in the mapping_file, specify "Input" for the sample in the "Strain" column. 
2. Specify the new library name in the "Library" column. 
3. Pass a new library_info.csv file containing the metadata for this new library. Crucially, the "Negative control cutoff" field in this library file can/should be empty (it will be ignored) for input files, because these libraries do not have their cutoffs determined yet.
4. Run the pipeline. Now, when running on any input samples, there will be a new file created, called input_samples.tsv. It has these fields:

```
FileName	Sample	Library	pGL0	Abundance	Cutoff
20230427-55-A10-MBS-SB	A10_MBS	pGL2_147	pGL0_12 [pAMβ1]	0.0015552099533437014	11.606399402107586
20230427-55-A10-MBS-SB	A10_MBS	pGL2_147	pGL0_113 [Dummy]	0.9984447900466563	23.158731726255038
20230427-55-A10-MBS-SB	A10_MBS	pGL2_147	pGL0_3 [RSF1010]	0.038880248833592534	12.038936025788614
20230427-55-A10-MBS-SB	A10_MBS	pGL2_147	pGL0_16 [2μ]	0.22706065318818042	14.219641503513795
20230427-55-A10-MBS-SB	A10_MBS	pGL2_147	pGL0_24 [pBBR1-UP]	1.401244167962675	27.826522790146132
```

Where the "Cutoff" column is the defined abundance cutoff. This is a z-test cutoff for the abundance (which is the reads of the origin normalized to dummy origin reads). It is calculated by assuming a standard deviation of 0.5 (the high range of observed standard deviations), a p-value of 0.05, and a bonferroni correction based on the number of origins in the library.

There are two test datasets of input libraries: `./test_data/test_inputs.csv` and `./test_data/test_inputs_and_strains.csv`, and their corresponding outputs. 
Note these aren't actual input libraries, they are just formatted as if they were (which is OK). Their samples are in amplicon_data.fastq.gz.

There is real output of our standard "v3" ORI libraries: `test_data/test_real_inputs.csv`. **This can be ingested and used for future conjugation runs of this library.**


## Whole plasmid sequencing ORI quantification

For identifying the presence of ORIs within a pool with whole plasmid (or whole genome) sequencing, use the script `whole_plasmid_quantification.py`. The inputs for this script are the directory of gzipped FASTQ files (a forward and reverse read file per sample) and a mapping file linking samples to their respective pools.

Example usage: `python whole_plasmid_quantification.py -d ./ms_data/fastq/ -m ./ms_data/mapping.csv`.

All arguments:

```
usage: whole_plasmid_quantification.py [-h] -d FASTQ_DIRECTORY -m MAPPING_FILE [-b BBMAP_FOLDER]

Whole plasmid sequencing of a plasmid ORI pool.

optional arguments:
  -h, --help            show this help message and exit
  -d FASTQ_DIRECTORY, --fastq_directory FASTQ_DIRECTORY
                        Directory of FASTQ files. File names must take the form: sample_*_R1_*.fastq.gz
  -m MAPPING_FILE, --mapping_file MAPPING_FILE
                        Mapping file of comma separated columns Sample,Pool.
  -b BBMAP_FOLDER, --bbmap_folder BBMAP_FOLDER
                        Directory containing BBTools on your system
```

Example outputs:

`./ms_data/plasmid_library_breadth.tsv` - Provides breadth of coverage for each ORI region in each sample. Recommend a breadth cutoff of at least 50% to call an ORI as present in a given sample.

`./ms_data/plasmid_library_coverage.tsv` - The mean coverage of each ORI region in each sample.
