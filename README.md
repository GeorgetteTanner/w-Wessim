

# w-Wessim
## Introduction

w-Wessim is an adapted version of Wessim (Kim S. et al., 2013), an *in silico* whole exome sequencing tool that combines the whole genome in silico sequencing tool, GemSIM (McElroy K. et al., 2012), with a BLAT alignment of exon capture hybridization probe sequences in order to define target exon regions for sequencing. w-Wessim builds upon Wessim to include modelling of copy number variants by weighting selection of exon capture probes based on the number of times they align to a genome, as well as allowing read lengths to be taken from a distribution to account for read trimming in the error model training data set. 

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

## Inputs

w-Wessim requires the following inputs:

* **Genome sequence:**
This is the genome that you intend to sequence.
FASTA format.

* **Probe sequences:**
These are the "probes" used in the BLAT alignment to define regions for w-Wessim to sequence. These can either be the sequences for exon capture kit hybridisation probes, or the sequences of real WES reads (recommended for more realistic read distributions from w-Wessim). 
Probe sequences for the Agilent SureSelect Human All Exon V4+UTRs kit and real WES reads (from the NCBI Sequence Read Archive, accession no. SRR2103613) that have been quality and adapter trimmed by cutadapt and filtered for a BWA MEM mapping quality (when aligned to hg38) of 60, are provided with w-Wessim. 

* **BLAT alignment of probe sequences:**
See the section below on BLAT.

* **GemSIM error model**
This tells w-Wessim how to incorporate errors in to the reads. Users can generate their own error models using GemSIM or use an existing model. One is provided with w-Wessim that has been trained on Illumina HiSeq 2000 WGS 101bp paired-end reads (from the NCBI Sequence Read Archive, accession no. ERR194146) from chromosome 1, that had been quality and adapter trimmed with cutadapt. Error models are specific to paired- or single-end read datasets and therefore the user must tell w-Wessim to generate paired end reads (with -p) if using a paired-end error model.


## BLAT

Blat is required to generate an alignment of the probe sequences to the genome the user wishes to sequence.

It takes a list of hybridisation probe sequences in FASTA format and a genome sequence, and outputs a .psl file listing locations of all alignments (that meet the given stringency parameters) for each probe.  

### Download

In the interest of speed, we reccomend using the multi-threaded pblat version of BLAT, created by Wang Meng and available from http://icebert.github.io/pblat/.
 
Alternatively, the single threaded original BLAT program can be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/ and selecting your computer type. 

fatotwobit is also required for use with (p)BLAT and can also be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/.

You may need to make the files executables by typing the following in the directory of the downloaded file: 

```
chmod a+x ./blat
chmod a+x ./fatotwobit
```


### Usage

It is often helpful to split the genome FASTA file that you wish to sequence into sections, or even individual chromosomes. This shortens the run time (if then ran in parallel) and also avoids errors that sometimes arrise when handling such large files. Splitting (and running in parallel), or subsampling, the probe sequences also shortens run time.

The genome (or sections) you want to sequence needs to be converted to 2bit format. This is done with:

```
faToTwoBit genome.fasta genome.fasta.2bit
```

(p)BLAT is then run with:

```
./blat genome.fasta.2bit probes.txt output.psl -minScore=95 -minIdentity=95
```

The -minScore and -minIdentity values can be varied to adjust how tight you want the alignment to be. We find that 95 for both gives the most realistic WES coverage results with w-Wessim and is also likely to accurately model the effect of variants on probe hybridisation in WES. See ```./blat```/```./pblat``` for more options.

### Combining .psl files

If the genome or probes were split when runnning (p)BLAT, the outputted .psl files need to be combined into one file and sorted for use with w-Wessim:

1.	Save the header:	
```
head -n 5 {any .psl file} > pslheader.txt
```

2. Remove the headers:
```
for f in *.psl ; do tail -n+6 $f > noheader$f ; done 
```

3. Combine noheader*.psl files: 
```
cat noheader*.psl > combined_noheader.psl
```

4. Sort combined file on column 10:
```
sort -k 10 -n combined_noheader.psl > sorted_combined_noheader.psl  
```

5. Add the header:
```
cat pslheader.txt  sorted_combined_noheader.psl > final.psl
```

## w-Wessim Usage

To run w-Wessim:

```
cd {w-Wessim_DIRECTORY}
python2 ./w-wessim2.py -R {INPUT_GENOME} -P {PROBES_FILE} -B {BLAT_OUTPUT}.psl -n 10000000 -l d -M {ERROR_MODEL}.gzip -pz -o {OUTPUT_NAME} -t 1 -v -m 20 -f 170 -d 35
```
### Parameters

|Parameter|Description|Default Value| 
|---|---|---|
|-R|Genome sequence in FASTA format - must be indexed with faidex. |Required
|-P|Probe sequences in FASTA format.|Required
|-B|Blat allignment output file of probes to the genome, in .psl format.|Required
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




