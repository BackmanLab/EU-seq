#!/bin/bash
##################################################################
# Author: Lucas Carter                                                   
# Email: lucascarter2025@u.northwestern.edu                                                          
# PI: Vadim Backman                                                   
# Description: 
# This script performs QC, trimming, mapping, and counting for
# EU-seq generated data. It generates directories, loops through 
# fastqs, and moves results into respective directories. This
# script is written to run on a SLURM HPC. A more thorough 
# description can be found at README.md
################################################################

##--------------------------------------------------------------## Module load

## Customize these for your environment
module load samtools/1.6
module load deeptools/3.1.1
module load fastqc
module load hisat2/2.1.0
module load TrimGalore/0.6.10

##--------------------------------------------------------------## Script begin

### Usage function tells users how to run the software
helpFunction()
{
      echo "*********************************** how to use NascentCov ***********************************"
      echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -t <threads> \n"
      echo -e "Run script in directory where fastqs are located \n"
      echo "This script completes 6 functions: Quality Control, Alignment, Sorting, Counting, and Read QC of features:"
      echo -e "\t -- Quality Control: eliminates the adaptor and low quality reads using trim_galore."
      echo -e "\t -- Mapping: mapping via HISAT2 for all features"
      echo -e "\t -- Sorting: Sort the bam files, keep uniquely mapping reads, generate coverage files  "
      echo -e "\t -- Counting: pile-up of reads on features in introns, protein coding, and non-coding RNAs \n \n"
      echo -e "\t -- Read QC: Check read quality of FASTQs at the end \n \n"
      echo -e "\t-h help \n"

      echo "For more detail information, please feel free to contact: lucascarter2025@u.northwestern.edu"
      echo "**************************************************"
   exit 1 # Exit script after printing help
}


##--------------------------------------------------------------## options

while getopts "f:w:g:r:t:" opt
do
   case $opt in
      f) R1=$OPTARG ;; ## forward read
      w) wdir=$OPTARG ;; ## working directory
      g) gendir=$OPTARG ;; ## path/to/bowtie_indices
      r) gtfdir=$OPTARG ;; ## path/to/gtf_file
      t) threads=$OPTARG ;; ## number of threads to HPC
      ?) echo "Unknown argument --${OPTARG}"; helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "${R1}" ] 
then
   echo "*** error: must specify forward read ***";
   helpFunction
fi

if [ -z "${wdir}" ] 
then
   echo "*** error: working directory where FASTQs are located must be provided ***";
   helpFunction
fi

if [ -z "${gendir}" ] 
then
   echo "*** error: path to reference genome indices must be provided ***";
   helpFunction
fi

if [ -z "${gtfdir}" ] 
then
   echo "*** error: path to genome annotation file must be provided ***";
   helpFunction
fi

if [ -z "${threads}" ] 
then
   echo "*** error: number of processors must be provided ***";
   helpFunction
fi

##--------------------------------------------------------------## Paths

# Get absolute path for scripts and results; check if required scripts exist

# Set working directory by moving script to directory where fastqs are located
cd ${wdir}

echo "$(date): Processing alignment of FASTQ sequence files"
echo -e "Generating directories. \n"

# make sure references are prepared correctly
root="$(dirname "$(pwd)")/"

echo "Workign directory is ${root}"
echo "Genome annotation file is located at ${gtfdir}"
echo -e "HISAT genome assembly is located in ${gendir} \n"

# Make .OUT directory for trimmed reads
outfold="${root}OUT"
[ ! -d $outfold ]&&mkdir $outfold

# make .SAM directory
samfold="${root}SAM"
[ ! -d $samfold ]&&mkdir $samfold

# Make QC directory
qcfold="${root}fastqc"
[ ! -d $qcfold ]&&mkdir $qcfold

# Make coverage directory
covfold="${root}coverage"
[ ! -d $covfold ]&&mkdir $covfold

# Make .BAM directory
bamfold="${root}BAM"
[ ! -d $bamfold ]&&mkdir $bamfold

# Make .HTSEQ_COUNTS directory
cntfold="${root}counts"
[ ! -d $cntfold ]&&mkdir $cntfold

##--------------------------------------------------------------## 1. Trim and QC
# begin data processing

cd ${wdir} ## move to fastq directory
R2=$(echo $R1 | sed 's/R1/R2/g')

###
echo "1. Starting the quality control: trimming ${R1} and ${R2}"
echo " **** 1.1 Running Fastqc on all FASTA files first: "

echo " **** 1.2 Running Trim Galore on all FASTA files now: "
echo -e " \t trim_galore --paired --retain_unpaired  --dont_gzip -o $outfold ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2}"

trim_galore --paired --retain_unpaired --dont_gzip -o $outfold ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2}

echo -e "${R1} and ${R2} have been processed \n"
###

##--------------------------------------------------------------## 2. Mapping

cd ${outfold} ## move to trimmed fastq folder
T1=$(echo $R1 | sed 's/R1_001.fastq.gz/R1_001_val_1.fq/g')
T2=$(echo $R2 | sed 's/R2_001.fastq.gz/R2_001_val_2.fq/g')
S1=$(echo $T1 | sed 's/R1_001_val_1.fq/R1_001.sam/g')

###
echo "2. Starting the mapping for ${T1} and ${T2}"
echo -e " \t hisat2 -p ${threads} -x ${gendir} -1 ${outfold}/${T1} -2 ${outfold}/${T2} -S ${samfold}/${S1} 2>${samfold}/hisat.err"
### Mapp reads using TopHat
hisat2 -p ${threads} -x ${gendir} -1 ${outfold}/${T1} -2 ${outfold}/${T2} -S ${samfold}/${S1} 2>${samfold}/hisat.err

echo -e "Mapping via HISAT2 complete for ${T1} and ${T2} \n"
###

##--------------------------------------------------------------## 3. Sort and statistics

cd ${samfold} ## move to SAM alignment folder
B1=$(echo $S1 | sed 's/R1_001.sam/filt.bam/g')
B2=$(echo $B1 | sed 's/filt.bam/filt.sort.bam/g')
B3=$(echo $B1 | sed 's/filt.bam/filt.n_sort.bam/g')

###
echo "3. sorting tophat mapping and retain the unique mapping reads for ${S1}"

samtools view -bS -F 1548 -q 30 ${samfold}/${S1} -o ${bamfold}/${B1} ## filter reads for unmapped and multimappers
samtools sort ${bamfold}/${B1} > ${bamfold}/${B2} ## Sort by position (for coverage file)
samtools sort -n  ${bamfold}/${B1} > ${bamfold}/${B3} ## Sort by read alignment flag
samtools index ${bamfold}/${B2} ## index reads
samtools flagstat ${bamfold}/${B2} > ${bamfold}/${B2}.FLAGSTAT.txt ## get alignment scores

echo -e "samtools finished processing ${S1} to ${B2} \n"
###

##--------------------------------------------------------------## 4. BAM coverage

BW1=$(echo $B2 | sed 's/_filt.sort.bam/.cpm.bigWig/g')
BW2=$(echo $B2 | sed 's/_filt.sort.bam/.rpgc.bigWig/g')

###
echo "4. Generating coverage files for ${B2}"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max"

# Bam coverage to generate BigWig
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max

echo -e "BAM coverage files generated for ${B2} using CPM and RPGC normalization \n"
### 

##--------------------------------------------------------------## 5. Counts

counts=$(echo $S1 | sed 's/_R1_001.sam/.counts/g')

###
echo "5. Counting reads for all annotations in ${B3}"
echo -e "\t htseq-count -f bam -q  -r name  -m union -s no  ${bamfold}/${B3} ${gtfdir} > ${cntfold}/${counts}"

htseq-count -f bam -q  -r name  -m union -s no  ${bamfold}/${B3} ${gtfdir} > ${cntfold}/${counts}

echo -e "Hiseq counting ended \n"

cd ${wdir} 
echo "6. QC on original fastq reads ${R1}"
echo -e "\t fastqc -t ${threads} ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2} -o ${qcfold}/"
ls -lh ${root}${PWD##*/}/${R1}
ls -lh ${root}${PWD##*/}/${R2}

fastqc -t ${threads} ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2} -o ${qcfold}/
echo -e " ******** All process have been completed ******** \n \n"
###




