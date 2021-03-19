#!/bin/bash

## call in the array number and pull out that line from the `lids` file, which should be a file with
## the base names of your libraries
LID=`sed -n ${SLURM_ARRAY_TASK_ID}p lids`

## trim R2s
module load contrib/trim_galore/0.4.5 ## For UW 
##module load trim_galore/0.4.0 ## For ASU ## But not sure about next line 
PATH=/sw/contrib/trim_galore/0.4.5:/sw/contrib/cutadapt/bin:$PATH

mkdir $PWD/trimmed_fastq

## Use if run as a script 
trim_galore $PWD/fastqs/${LID}.R2.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $PWD/trimmed_fastq

### Worked for 1 at UW 
###trim_galore $PWD/fastqs/LID_104229.R2.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $PWD/trimmed_fastq


##align reads using kallisto
module load contrib/kallisto/0.46.1
mkdir $PWD/kallisto_alignment

kallisto_index=/gscratch/csde/smacklab/genomes/kallisto/Mmul_10/mmul10
input_folder=$PWD/trimmed_fastq
kallisto_output_folder=$PWD/kallisto_alignment

kallisto quant -i $kallisto_index -t 24 -o $kallisto_output_folder/${LID} --single-overhang --single -l 200 -s 200 $input_folder/${LID}.R2_trimmed.fq.gz 

### Worked for 1 
#kallisto quant -i $kallisto_index -t 24 -o $kallisto_output_folder/LID_104229 --single-overhang --single -l 200 -s 200 $input_folder/LID_104229.R2_trimmed.fq.gz


## calculate stats
fastq_folder=$PWD/fastqs/
trimmed_fastq=$PWD/trimmed_fastq/

echo $LID, $(expr $(zcat $fastq_folder/$LID*R2*.fastq.gz|wc -l) / 4), $(expr $(zcat $trimmed_fastq/$LID*R2*.gz|wc -l) / 4), $(awk '{s+=$4}END{print s}' $kallisto_output_folder/$LID/abundance.tsv), $(sed -n '7p' $kallisto_output_folder/$LID/run_info.json | sed 's/.*://g' | sed 's/,//g')

