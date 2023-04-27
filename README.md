sciCTextract
============

Simple sciCUT&amp;Tag demultiplexing.

Installation
------------

Installation within a conda environment or virtualenv is recomended.
In an active enviroment with python3 installed, clone the git
repository and install with pip:
```
pip install .
```
This will install dependencies such as Biopython and put the
`sciCTextract` command on your path.

Input Files
-----------

Demultiplexing requires four input Fastq files following Illumina naming
conventions, including read headers of the form:
```
@VH00319:342:AACKYJMM5:1:1101:31410:1000 1:N:0:0
```
On multi-lane flowcells, reads should not be split by lane, unless needed
to process lanes independently (to support XP worflows for example). A
demultiplexing run requires exactly four Fastq files for paired sequence
reads (_R1 &amp; _R2) and paired index reads (_I1 &amp; _I2).

The process also requires two barcode tables in comma-sepearated value
format, one for Tn5 barcodes that will define the prefixes of the samples
names and one for the Primer barcodes that will defined the suffixes.

The Tn5 barcode table should have this form at minimum:

|  Sample Name  | Tn5\_s7   | Tn5\_s7\_seq | Tn5\_s5   | Tn5\_s5\_seq |
| :-----------: | :-------: | :----------: | :-------: | :----------: |
| Hs\_H3K27ac   | P7\_i7\_1 |  ATTACTCG    | P5\_i5\_1 |  TATAGCCT    |
| Hs\_H3K27ac   | P7\_i7\_1 |  ATTACTCG    | P5\_i5\_2 |  ATAGAGGC    |

The Primer barcode table should have this form at minimum:

| i7\_index\_seq | i5\_index\_seq | i7\_index\_id | i5\_index\_id |  ID  |
| :------------: | :------------: | :-----------: | :-----------: | :--: |
|   GGACTCCT     |    TAGATCGC    |      P7\_5    |      P5\_1    | 10pM |
|   GGACTCCT     |    CTCTCTAT    |      P7\_5    |      P5\_2    | 10pM |

Note that barcodes should _always_ be specified in forward strand
("Workflow A") orientation. This allows the same barcode tables to be
used with different types of Illumina instruments.
All barcodes are currently required to be 8nt.

Running
--------

With four Fastq files in hand and two barcode tables defined,
create an output directory (e.g., `mkdir fastq_out`) and launch
demultiplexing, for example:

```
sciCTextract \
    --outdir fastq_out \
    --Tn5_Barcode Tn5_Barcode_Annotation.csv \
    --Primer_Barcode Primer_Barcode_Annotation.csv \
    Undetermined_S0_R1_001.fastq.gz \
    Undetermined_S0_R2_001.fastq.gz \
    Undetermined_S0_I1_001.fastq.gz \
    Undetermined_S0_I2_001.fastq.gz
```

Note that our current typical use is to run on the Illumina NextSeq 2000.
Default settings should work for NextSeq 1000/2000, NovaSeq 6000 (v1.5 or more
recent). For instruments that use forward-strand workflows
(MiSeq, HiSeq, MiniSeq Rapid, etc.) we provide the `--forward-mode` option
to override the default reverse-complementing of the i5 barcode reads.

Output
------

Output consists of one pair of gzip compressed Fastq files per sample.
The read headers are re-written to include the error-corrected barcode
sequences and to be compact while retaining enough information to
unambihguously identify each source read. For example:

```
@HMH53BCX3:1:1105:11433:2512\_GCGTTAAA\_GTGTATCG\_AGCGATAG\_CAGGACGT 1:N:0:0
```
