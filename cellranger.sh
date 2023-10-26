#!/usr/bin/bash

# Job Name
#$ -N cellranger_multi

# Execute the script from the Current Working Directory
#$ -cwd

# Make sure your enviornment variables are preserved
#$ -V

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'sge-output' in the current working directory (cwd)
#$ -o cellranger_multi_all

#Specify number of threads
#$-pe threaded 8

# Tell the job your memory requirements
#$ -l h_vmem=15G 

# Send mail when the job is submitted, and when the job completes (n means no)
#$ -m n

#  Specify an email address to use
#$ -M foo.cheung@nih.gov

#Array number of jobs (count number of files with ls | wc -l)
#$ -t 12 

# Specify location of the files to be processed


# Create a bash array of numbers for the files 
SAMPLE_LIST_R1=(1 2 3 4 5 6 7 8)


# Get index from $SGE_TASK_ID
INDEX=$((SGE_TASK_ID-1))

module load cellranger/7.1.0
cellranger multi --id=NEW4_multi_config_${SAMPLE_LIST_R1[$INDEX]} --csv=multi_config_${SAMPLE_LIST_R1[$INDEX]}.csv  --localcores=8 --localmem=120

