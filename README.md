## w-Wessim2 is now available at https://github.com/GeorgetteTanner/w-Wessim2
Improvements to w-Wessim2 include using ReSeq (Schmeing S and Robinson M, 2020) for improved sequencing error modelling, as well as requiring less computational resources than w-Wessim.

# w-Wessim

## Introduction

w-Wessim is an adapted version of Wessim (Kim S. et al., 2013), an *in silico* whole exome sequencing (WES) tool that combines a method for selecting fragments from target regions with the whole genome in silico sequencing tool, GemSIM (McElroy K. et al., 2012). Target regions are determined through a BLAT alignment of exon capture hybridization probe sequences to the genome to be sequenced. w-Wessim builds upon Wessim to include modelling of copy number variants by weighting selection of exon capture probes based on the number of times they align to a genome, as well as allowing read lengths to be taken from a distribution to account for read trimming in the error model training data set. 

We also provide an improved protocol that more accurately models the read distributions seen in real WES data by using real WES read sequences instead of exon capture probe sequences for the BLAT alignment.

NB. The Error model provided with w-Wessim (trained using GemSIM), has a higher than expected base quality score distribution, specifically for error bases. This may cause a high false positive rate during point variant calling, which users can overcome either by training their own error models with GemSIM, or by subsampling false positive calls to more realistic proportions. 

For more information on w-Wessim see: https://doi.org/10.1093/bioinformatics/bty1063

Please cite the following when using w-Wessim:

Tanner G, Westhead DR, Droop A, Stead LF. Simulation of heterogeneous tumour genomes with HeteroGenesis and in silico whole exome sequencing. Bioinformatics. 2019 Jan 4;35(16):2850–2. 

...along with the authors of Wessim and GemSIM:

Kim S, Jeong K, Bafna V. Wessim: a whole-exome sequencing simulator based on in silico exome capture. Bioinformatics. 2013;29:1076–7.

McElroy KE, Luciani F, Thomas T. GemSIM: general, error-model based simulator of next-generation sequencing data. BMC genomics. 2012;13(1):74.

For any issues in running w-wessim, please leave a comment on github or email medgnt@leeds.ac.uk.

## Versions

v2 - Improved resource management and speed. (25/06/2020)

v1 - Initial commit. (09/08/2018)

## Requirements

w-Wessim requires Python2, with numpy and pysam (which requires htslib and samtools). w-Wessim has been tested with the following:

Python 2.7.12

numpy 1.12.0

w-Wessim can either be run with probe sequences from exon capture kits, or with real WES reads as probes. A BLAT alignment of either of these probe sets to the genome being sequenced is necessary prior to running w-Wessim.

The first option, with exon capture kit probes, requires relatively low computational requirements and can be done on a standard computer within a few hours.

The second option, with real WES reads as probes, produces more realistic coverage across target and off target regions but requires a high performance computing system due to both time and memory requirements. The real read probes should be randomly downsampled, depending on the number of simulated reads needing to be generated; Using at least the same number of probes as the number of reads being simulated allows for the most realistic distribution overall, but requires high resources. Using fewer probe numbers (eg. ten times fewer probes than the number of reads needing to be generated) reduces resource requirements, but may result in clumping of the off-target reads, although it has less impact on on-target regions where there's a higher density of real WES probes. 

Using the full set of 1x10^8 real read probes used in the example in the paper, requires 5GB RAM and ~200h for the BLAT alignment on a single thread (this can be multithreaded, or ran across multiple nodes in separate runs by splitting the read number, probes or genome sequence to shorten run time) and around 50GB RAM for w-Wessim, to sequence the hg38 human reference genome. Whereas, using 1x10^7 requires around 20h and 5GB RAM for the BLAT on a single thread, and 7GB RAM by w-Wessim. The time required for w-Wessim mostly depends on the number of reads being generated; 1x10^7 reads requires around 8h. Larger genomes (eg. when using diploid human genomes) will require higher resources.

## Installation
```bash
git clone https://github.com/GeorgetteTanner/w-Wessim.git
```

## Inputs

w-Wessim requires the following inputs:

* **Genome sequence:**
This is the genome that you intend to sequence.
FASTA format. Must be indexed with faidx.

* **Probe sequences:**
These are the "probes" used in the BLAT alignment to define regions for w-Wessim to sequence. These can either be the sequences for exon capture kit hybridisation probes, or the sequences of real WES reads (recommended for more realistic read distributions from w-Wessim). Probe names must be converted into integers prior to running the BLAT. Probe sequences for the Agilent SureSelect Human All Exon V4+UTRs kit (or any other kit for which probe sequences are avaialable) can be downloaded from https://earray.chem.agilent.com/suredesign/index.htm and converted to FASTA format with the Prep\_Probe2Fa.py script from the orginal Wessim tool (http://sak042.github.io/Wessim/). Real WES reads (from the NCBI Sequence Read Archive, accession no. SRR2103613, captured with the Agilent SureSelect Human All Exon V5+UTRs kit) that have been quality and adapter trimmed by cutadapt and filtered for a high BWA MEM mapping quality, are provided as sample ERR2752113 from the European Nucleotide Archive (http://ftp.sra.ebi.ac.uk/vol1/run/ERR275/ERR2752113/real_wes_reads_probes.fastq.gz). These will need to be converted from fastq to fasta format with:

```
paste - - - - < file.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa
```
 These can be downsampled to the required number (eg. 10000000) with:

```
cat file.fa| paste - - > shuf -n 10000000 | sed 's/^.//' > shuf_file.fa
sort -k1 -n shuf_file.fa > sort_file.fa 
sed -e 's/^/>/' sort_file.fa | tr '\t' '\n' > downsampled_10m_file.fa
```



* **BLAT alignment of probe sequences:**
See the section below on BLAT.

* **GemSIM error model**
This tells w-Wessim how to incorporate errors in to the reads. Users can generate their own error models using GemSIM or use an existing model. One is provided with w-Wessim ('lib/hs2000p.gzip') that has been trained on Illumina HiSeq 2000 WGS 101bp paired-end reads (from the NCBI Sequence Read Archive, accession no. ERR194146, and quality and adapter trimmed with cutadapt) from chromosome 1, with known variant sites ignored. 

If using an error model, such as the one provided, that has been trained on cleaned data it may not be appropriate to clean the resulting reads from w-Wessim.

## BLAT

BLAT is required to generate an alignment of the probe sequences to the genome the user wishes to sequence.

It takes a list of hybridisation probe sequences in FASTA format and a genome sequence, and outputs a .psl file listing locations of all alignments (that meet the given stringency parameters) for each probe.  

### Download

In the interest of speed, we reccomend using the multi-threaded pblat version of BLAT, created by Wang Meng and available from http://icebert.github.io/pblat/.

Alternatively, the single threaded original BLAT program can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/ and selecting your computer type. 

faToTwoBit is also required for use with (p)BLAT and can also be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/.

You may need to make the files executables by typing the following in the directory of the downloaded file: 

```
chmod a+x ./blat
chmod a+x ./fatotwobit
```


### Usage

Probe names must be converted into integers prior to running the BLAT in order to work with w-Wessim. This is already done for the available real WES read probes.

The genome you want to sequence needs to be converted to 2bit format. This is done with:

```
faToTwoBit genome.fasta genome.fasta.2bit
```

(p)BLAT is then run with:

```
blat genome.fasta.2bit probes.txt output.psl -minScore=95 -minIdentity=95
#or
pblat genome.fasta.2bit probes.txt output.psl -threads=8 -minScore=95 -minIdentity=95

```

The -minScore and -minIdentity values can be varied to adjust how tight you want the alignment to be. We find that 95 for both gives the most realistic WES coverage results with w-Wessim and is also likely to accurately model the effect of variants on probe hybridisation in WES. See ```./blat```/```./pblat``` for more options.


## w-Wessim Usage

To run w-Wessim:

```
cd {w-Wessim_DIRECTORY}
python2 w-wessim.py -R {INPUT_GENOME} -P {PROBES_FILE} -B {BLAT_OUTPUT}.psl -n 10000000 -l d -M {ERROR_MODEL}.gzip -o {OUTPUT_NAME} -m 20 -f 170 -d 35
```
### Parameters

|Parameter|Description|Default Value| 
|---|---|---|
|-R|Genome sequence in FASTA format - must be indexed with faidex. |Required
|-P|Probe sequences in FASTA format.|Required
|-B|BLAT allignment output file of probes to the genome, in .psl format.|Required
|-M|GemSIM error model.|Required
|-o|Output file name - .fasta(.gz) will be attached to this.|Required
|-n|Number of reads/pairs of reads.|Required
|-l|Read length. 'd' may be given instead of an integer for read lengths to be taken from the distribution provided in the GemSIM error model.|Required
|-z|Compress output with gzip.|[false]
|-f|Mean fragement size (when using paired-end sequencing).|200
|-d|Standard deviation of fragment size.|50
|-m|Minimum fragment lenngth. |read length + 20 
|-y|Minimum required fraction of probe match to be hybridized.|50
|-w|Penalty weight for indel in the hybridization.|2
|-q|Quality score offset.|33



## Example


This example demonstrates how to use w-Wessim with the real reads probe set to create the most realistic sequencing dataset. This is a very time and memory consuming process and may not be feasible without access to a high performance computing system. Instead, the option of using a subsampled probe set (with 1/1000th of the probes) is available if the user wishes to run the programs on a standard computer just for testing the code. The resulting sequencing data set from this will look very patchy and is not intended for use.
Alternatively, the probe sequences from an exon capture kit can be used. This is much quicker and less memory intensive but results in a slightly less realistic distribution of reads. 

Instructions for all three options are included below.


```bash
####Download programs:

git clone https://github.com/GeorgetteTanner/w-Wessim.git

#(need to find the correct binary file for your operating system:)
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit
chmod a+x ./fatotwobit

#1.
git clone https://github.com/icebert/pblat.git
cd pblat
make
cd ..
#OR 2. (if conda is installed)
conda install pblat

####Download probes:

#EITHER:

#1. Download full set of real WES reads probe sequences and convert from fastq to fasta:
cd w-Wessim
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA157/ERA1574375/fastq/real_wes_reads_probes.fastq.gz
gunzip real_wes_reads_probes.fastq.gz
paste - - - - < real_wes_reads_probes.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > real_wes_reads_probes.fa

#OR

#2.Download 1/1000th of real WES reads probe sequences:
cd w-Wessim
wget https://github.com/GeorgetteTanner/data/raw/master/real_wes_reads_probes_subsampled.fa.gz
gunzip real_wes_reads_probes_subsampled.fa.gz
mv real_wes_reads_probes_subsampled.fa real_wes_reads_probes.fa

#OR

#3. Manually download the Agilent SureSelect Human All Exon V4+UTRs 
#kit probes from https://earray.chem.agilent.com/suredesign/index.htm 
#and convert to FASTA format with the Prep\_Probe2Fa.py script from 
#the orginal Wessim tool (http://sak042.github.io/Wessim/). 

#code not shown here for this#


####Download a reference genome or any other genome you wish to sequence:
mkdir ../test1
cd ../test1
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
gunzip Homo_sapiens_assembly38.fasta.gz

###BLAT
#(Older versions of pblat had issues with large genomes. If this happens, you can split the genome
#and run pblat on each section separately. This doesn't appear to be a problem with newer pblat versions.)

faToTwoBit Homo_sapiens_assembly38.fasta Homo_sapiens_assembly38.fasta.2bit
pblat Homo_sapiens_assembly38.fasta.2bit ../w-Wessim/real_wes_reads_probes.fa blatoutput_Homo_sapiens_assembly38.psl -minScore=95 -minIdentity=95 -threads=8

###w-Wessim:
cd ../w-Wessim
python2 w-wessim.py -R ../test1/Homo_sapiens_assembly38.fasta -P real_wes_reads_probes.fa -B ../test1/blatoutput_Homo_sapiens_assembly38.psl -n 100000 -l d -M lib/hs2000p.gzip -o ../test1/w-wessimoutput -m 20 -f 170 -d 35

```
