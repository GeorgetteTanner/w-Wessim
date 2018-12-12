

# w-Wessim
## Introduction

w-Wessim is an adapted version of Wessim (Kim S. et al., 2013), an *in silico* whole exome sequencing (WES) tool that combines a method for selecting fragments from target regions with the whole genome in silico sequencing tool, GemSIM (McElroy K. et al., 2012). Target regions are determined through a BLAT alignment of exon capture hybridization probe sequences to the genome to be sequenced. w-Wessim builds upon Wessim to include modelling of copy number variants by weighting selection of exon capture probes based on the number of times they align to a genome, as well as allowing read lengths to be taken from a distribution to account for read trimming in the error model training data set. 

We also provide an improved protocol that more accurately models the read distributions seen in real WES data by using real WES read sequences instead of exon capture probe sequences for the BLAT alignment.

For more information, see "Simulation of Heterogeneous Tumour Genomes with HeteroGenesis and In Silico Whole Exome Sequencing, Tanner G et al., 2018." (Manuscript submitted)
Please cite this when using w-Wessim, along with the authors of Wessim and GemSIM:

Kim S, Jeong K, Bafna V. Wessim: a whole-exome sequencing simulator based on in silico exome capture. Bioinformatics. 2013;29:1076–7.

McElroy KE, Luciani F, Thomas T. GemSIM: general, error-model based simulator of next-generation sequencing data. BMC genomics. 2012;13(1):74.


## Requirements

w-Wessim requires Python2, with numpy and pysam (which requires htslib and samtools). w-Wessim has been tested with the following:

Python 2.7.12

numpy 1.12.0

pysam 0.10.0

samtools 1.3.1

htslib 1.3.2

w-Wessim can either be run with probe sequences from exon capture kits, or with real WES reads as probes. 

The first option requires relatively low computational requirements and can be done on a standard computer within a few hours.

The second option produces more realistic coverage across target and off target regions but requires a high performance computing system due to both time and memory requirements - around 9h/threads and 72GB RAM x threads to generate 1 × 10^7 pairs of reads. The large memory requirement results from the high number of real reads used as probes, which can be downsampled if needed. A BLAT alignment of the probes to the genome being sequenced is necessary prior to running w-Wessim. For the full set of 1x108 real read probes used in our example, this would take ~700h/threads. However, this can be ran across multiple nodes in separate runs by splitting the read number, probes or genome sequence.

##Installation
```bash
git clone https://github.com/GeorgetteTanner/w-Wessim.git
```

## Inputs

w-Wessim requires the following inputs:

* **Genome sequence:**
This is the genome that you intend to sequence.
FASTA format. Must be indexed with faidx.

* **Probe sequences:**
These are the "probes" used in the BLAT alignment to define regions for w-Wessim to sequence. These can either be the sequences for exon capture kit hybridisation probes, or the sequences of real WES reads (recommended for more realistic read distributions from w-Wessim). Probe sequences for the Agilent SureSelect Human All Exon V4+UTRs kit (or any other kit for which probe sequences are avaialable) can be downloaded from https://earray.chem.agilent.com/suredesign/index.htm and converted to FASTA format with the Prep\_Probe2Fa.py script from the orginal Wessim tool (http://sak042.github.io/Wessim/). Real WES reads (from the NCBI Sequence Read Archive, accession no. SRR2103613, captured with the Agilent SureSelect Human All Exon V5+UTRs kit) that have been quality and adapter trimmed by cutadapt and filtered for a high BWA MEM mapping quality, are provided as sample ERR2752113 from the European Nucleotide Archive (ftp://ftp.sra.ebi.ac.uk/vol1/ERA157/ERA1574375/fastq/real\_wes\_reads\_probes.fastq.gz). These will need to be converted from fastq to fasta format with:

```
paste - - - - < file.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > file.fa
```

* **BLAT alignment of probe sequences:**
See the section below on BLAT.

* **GemSIM error model**
This tells w-Wessim how to incorporate errors in to the reads. Users can generate their own error models using GemSIM or use an existing model. One is provided with w-Wessim ('lib/hs2000p.gzip') that has been trained on Illumina HiSeq 2000 WGS 101bp paired-end reads (from the NCBI Sequence Read Archive, accession no. ERR194146, and quality and adapter trimmed with cutadapt) from chromosome 1, with known variant sites ignored. Error models are specific to paired- or single-end read datasets and therefore the user must tell w-Wessim to generate paired end reads (with -p) if using a paired-end error model.

If using an error model, such as the one provided, that has been trained on cleaned data it is not recommended to clean the resulting reads from w-Wessim.

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

It's often helpful to split the genome FASTA file that you wish to sequence into sections to avoid errors that sometimes arrise with pblat when handling larger files - we find splitting a human genome into 4 files (each ~1.5GB) is sufficient. However, each run still takes a similar length of time regardless of how large the genome section is, so its best to not split the genome up too much. Splitting the probe sequences and running them in parallel on multiple nodes is also possible if wanting to shorten run time.

The genome (or sections) you want to sequence needs to be converted to 2bit format. This is done with:

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
python2 w-wessim2.py -R {INPUT_GENOME} -P {PROBES_FILE} -B {BLAT_OUTPUT}.psl -n 10000000 -l d -M {ERROR_MODEL}.gzip -pz -o {OUTPUT_NAME} -t 1 -v -m 20 -f 170 -d 35
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
|-p|Generate paired-end reads - must match the GemSIM error model.|[false]
|-z|Compress output with gzip.|[false]
|-f|Mean fragement size (when using paired-end sequencing).|200
|-d|Standard deviation of fragment size.|50
|-m|Minimum fragment lenngth. |read length + 20 
|-w|Minimum required fraction of probe match to be hybridized.|50
|-y|Penalty weight for indel in the hybridization.|2
|-t|Number of threads.|1
|-q|Quality score offset.|33



## Example


This example demonstrates how to use w-Wessim with the real reads probe set to create the most realistic sequencing dataset. This is a very time and memory consuming process and not feasible without access to a high performance computing system. The option of using a subsampled probe set (with 1/1000th of the probes) is available if the user wishes to run the programs on a standard computer just for testing the code. The resulting sequencing data set from this will look very patchy and is not intended for use.
Alternatively, the probe sequences from an exon capture kit can be used. This is much quicker and less memory intensive but results in a slightly less realistic distribution of reads. 

Instructions for all three options are included below.


```bash
#Download programs - (you may get a few warnings during pblat installation that can be ignored):
git clone https://github.com/GeorgetteTanner/w-Wessim.git
git clone https://github.com/icebert/pblat.git
#(need to find the correct binary file for your operating system:)
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToTwoBit
chmod a+x ./fatotwobit

#EITHER:

#1. Download full set of probe sequences and convert from fastq to fasta:
cd w-Wessim
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA157/ERA1574375/fastq/real_wes_reads_probes.fastq.gz
gunzip real_wes_reads_probes.fastq.gz
paste - - - - < real_wes_reads_probes.fastq | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > real_wes_reads_probes.fa

#OR

#2.Download subsampled probe sequences:
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


#Download a reference genome or any other genome you wish to sequence:
mkdir ../test1
cd ../test1
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai
gunzip Homo_sapiens_assembly38.fasta.gz

#Split genome into chromosomes (ignoring unplaced scafolds):
for chr in $(seq 1 22) X Y; do samtools faidx Homo_sapiens_assembly38.fasta ${chr} > Homo_sapiens_assembly38_${chr}.fasta ; done

#Combine fasta files. These need to be grouped into files of no more 
#than around 2GB to avoid errors with pblat - we recommend grouping 
#into 2 for haploid genomes or 4 for diploid genomes. Each run of
#pblat takes a similar length of time regardless of #genome length so 
#its best to group into as few files as possible.: 

#EITHER:

#1. For haploid:
cat Homo_sapiens_assembly38_chr1*fasta > Homo_sapiens_assembly38_1.fasta
cat Homo_sapiens_assembly38_chr[^1]*fasta > Homo_sapiens_assembly38_2.fasta

#OR

#2. For diploid (which requires further splitting to reduce file sizes):
cat Homo_sapiens_assembly38_chr1*A*fasta > Homo_sapiens_assembly38_1.fasta
cat Homo_sapiens_assembly38_chr1*B*fasta > Homo_sapiens_assembly38_2.fasta
cat Homo_sapiens_assembly38_chr[^1]*A*fasta > Homo_sapiens_assembly38_3.fasta
cat Homo_sapiens_assembly38_chr[^1]*B*fasta > Homo_sapiens_assembly38_4.fasta 

#Convert each chromosome fasta to 2bit:
for f in Homo_sapiens_assembly38_*.fasta ; do faToTwoBit $f ${f}.2bit ; done

#pblat (change probe directory and name if required): 
for f in Homo_sapiens_assembly38_*.fasta.2bit ; do ../pblat/pblat $f ../w-Wessim/real_wes_reads_probes.fa blatoutput_$(basename $f .fasta.2bit).psl -minScore=95 -minIdentity=95 -threads=8; done

#Combine .psl files:
#Save the header:
head -n 5 $(ls prefix*.fasta | head -1 ) > pslheader.txt
#Remove the headers:
for f in *.psl ; do tail -n+6 $f > noheader_$f ; done 
#Combine noheader*.psl files and sort combined file on column 10: 
for clone in clone1 clone2 germline ; do cat noheader_${clone}*.psl | sort -k 10 -n > sorted_combined_noheader_${clone}.psl ; done
#Add the header:
cat pslheader.txt sorted_combined_noheader_${clone}.psl > ${clone}.psl

#Combine fasta files into full genomes:
for clone in clone1 clone2 germline ; do cat prefix${clone}*.fasta > prefix${clone}.fasta ; done

#w-Wessim:
cd ../w-Wessim
python2 w-wessim2.py -R ../test1/prefix${clone}.fasta -P real_wes_reads_probes.fa -B ${clone}.psl -n 100000 -l d -M lib/hs2000p.gzip -pz -o ../test1/w-wessimoutput -t 1 -v -m 20 -f 170 -d 35

```