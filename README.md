# CRISPRone - command line version

> *Warning: The command-line version of this software tool was not originally developed for public consumption, and **output formats may not be in a friendly human-readable format**. Command-line version is provided here to allow interested parties who may find interest in utilizing the command-line version rather than relying on the [web-based server](https://omics.informatics.indiana.edu/CRISPRone/).*

CRISPRone provides annotation of CRISPR—Cas systems including:

- CRISPR arrays of repeat-spacer units, and cas genes
- type (and subtype) of predicted system(s)
- anti-repeats (part of tracrRNA genes in type II CRISPR–Cas systems)

## Prerequisites & Dependencies
Dependencies required to be installed prior to running CRISPRone.
```
    - Python 2.7+ OR Python 3.6+
    - Perl 
```

## Usage
```
usage: crisprone-local.py [output-base] [output-prefix] [input-fna] <input-gff> 

required arguments:
  output-base   Output directory for CRISPRone results.
  output-prefix Prefix/label for output subdirectories and output files.
  intput-fna    Path to input FASTA file.

optional arguments:
  input-gff     Path to GFF file of intput FASTA file.
```

## Citation
[Quan Zhang and Yuzhen Ye. Not all predicted CRISPR-Cas systems are equal: Isolated cas genes and classes of CRISPR like elements. BMC Bioinformatics, 18:92, 2017](https://www.readcube.com/articles/10.1186%2Fs12859-017-1512-4?author_access_token=Qirj1Fc54XCKXl2D5HFJCm_BpE1tBhCbnbw3BuzI2RO4ZMi96tlS0oMzwlg-pp16kdzldidFLmNumT7rd8h_qbSZ5oHFl3YSgoHGATPEpRriHC-flZ89ve1ZvCelk1nsA5g-3ePZ_RFmWnd1b8Tjyw==)

## Acknowledgement
This work was supported by NSF grants DBI-0845685 & DBI-1262588
