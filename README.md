# bellerophon

*——named after the Greek hero who slayed the chimera*

### Introduction

This tool was made for analysing MHC typing and other low-complexity gene amplicon data (as opposed to high-complexity data such as whole-microbiome 16S rRNA amplicons) generated on PacBio system. It performs allele calling while detecting polymorphic sites within the sequences and removing potential chimeric sequence variants.

The data processing method shown in the flowchart below was developed for high-resolution allele calling, taking into account that technical errors and artefacts such as PCR heteroduplexes and chimeras can be severe in mixed-template or multi-target amplifications (common in MHC typing). Briefly, circular consensus sequence (CCS) reads are generated, demultiplexed, trimmed, filtered, and aligned using tools available from PacBio software suite SMRTLink. With processed reads aligned to a reference sequence, bellerophon first generates a multiple sequence alignment with the reads, then identifies real polymorphic sites within the sequences, and then calls real alleles and putative chimeric sequences based on read support. Further details are available from [*(placeholder: link to paper)*](http://). 

![flowchart](img/flowchart.png)

### How to use

Simply clone this repository into a folder where you want to save it:
```
git clone https://github.com/yuanyuan929/bellerophon.git
```
You will also need Python3, Biopython, and samtools.

Try running bellerophon using the test data provided:
```
cd bellerophon/toy_data
python ../bellerophon.py sample.bam ref.fa --db db.fa --blacklist blacklist.txt
```
Example output files can be found in bellerophon/toy_data/output_files


### Documentation

***bellerophon.py***

This is the allele caller, the main tool of this repository.

```
usage: bellerophon.py [-h] [--db known_allele_db] [--blacklist blacklist] [--min_var_freq <float>] [--min_read_perc <float>] [--evidence] input_bam reference_seq

Call alleles from gene amplicon data generated on PacBio system.

positional arguments:
  input_bam             bam file produced by pbalign (required)
  reference_seq         fasta file containing the reference sequence used for pbalign (required)

optional arguments:
  -h, --help            show this help message and exit
  --db known_allele_db  fasta file containing known alleles (optional)
  --blacklist blacklist
                        region(s) in the sequences to exclude from variable site calling (optional)
  --min_var_freq <float>
                        minimum frequency required to call nucleotide variant (default: 0.05)
  --min_read_perc <float>
                        minimum percentage of reads required to consider a sequence variant as a candidate (default: 0.01)
  --evidence            generates an additional output file listing read IDs of each allele
```
`--db` If there are known alleles available for the genes of interest, the sequences of known alleles can be provided in a fasta file. Bellerophon will compare identified alleles with the known alleles and assign allele names if matches are found.

`--blacklist` If the target gene contains regions comprising long streches of a single type of nucleotides, these homopolymeric regions may need to be blocked due to high rates of indel errors caused during PCR & sequencing which can interfere with allele calling. Usually these segments are more likely to occur in noncoding regions than in coding regions ([ref1](http://dirac.cnrs-orleans.fr/~piazza/PB/files/DNA.pdf), [ref2](https://academic.oup.com/nar/article/26/17/4056/1176756)). To use this option, provide the start and end coordinates (1-based; corresponding to positions in the reference sequence) for regions to be blocked in a tab-delimited file, with each region on a separate line ([example](toy_data/blacklist.txt)).

`--min_var_freq` This number is used for distinguishing between sequencing errors and real variations. A nucleotide variant with a frequency among all reads that is higher than this cut-off value is considered a real variant. This value should be adjusted based on the number of genes co-amplified: the more genes co-amplified, the lower the value (shouldn't be lower than the error rate in the reads though).

`--min_read_perc` This is for filtering sequence variants based on relative read abundance. E.g. 0.05 means a sequence variant must be supported by at least 5% of the reads to be considered a candidate.

`--evidence` Use this option if you would like to know the read IDs for each allele called.


***bam_to_aligned_fasta.py***

This tool converts aligned reads in bam format into a multiple sequence alignment in fasta format.

```
usage: bam_to_aligned_fasta.py [-h] [--min_insert_freq <float>] [--noInsertionFilter] input_bam reference_seq

Generate multiple sequence alignments in fasta format from bam file.

positional arguments:
  input_bam             bam file produced by pbalign (required)
  reference_seq         fasta file containing the reference sequence used for pbalign (required)

optional arguments:
  -h, --help            show this help message and exit
  --min_insert_freq <float>
                        minimum frequency required to call real insertions (default: 0.2)
  --noInsertionFilter   keep all insertions (default is to remove putatively artefactual insertions based on read frequency)
```
`--min_insert_freq` The threshold value for filtering real insertions from indel errors. Due to relatively high indel rates in PacBio reads, a higher value is recommended.

`--noInsertionFilter` All insertions will be kept (the resulting multiple sequence alignment will likely contain lots of gaps).

***filter_fasta_by_sequence_id.py***

This tool is for extracting or removing sequences from a fasta file. Input file `Sequence_ID_list` should be in the format of one ID per line.

```
usage: filter_fasta_by_sequence_id.py [-h] [--remove] Input_fasta Sequence_ID_list Output_fasta

Extract or remove sequences from a fasta file by sequence IDs.

positional arguments:
  Input_fasta       Input fasta file (required)
  Sequence_ID_list  Input file containing sequence IDs (required)
  Output_fasta      Name of output fasta file (required)

optional arguments:
  -h, --help        show this help message and exit
  --remove          Remove listed sequences (Default: extract listed sequences)
```

### Citation

Please cite: [*(placeholder: link to paper)*](http://)

