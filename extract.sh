#!/bin/bash

# Specifications for eddie

#$ -N Single_read_extract                     # Name the job
#$ -pe sharedmem 32                              # Request 32 core cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=12G                                 # Request 12G memory per cpu
#$ -cwd
#$ -l h_rt=48:00:00                              # Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh
export PATH="/gpfs/igmmfs01/eddie/sproul-lab/Gawain/tools:$PATH" #load the modkit tool
# Arguments:
# $1: path to the modifiedbam file
# $2: output path for tsv contain single read methylation information
# $3: path to reference genome
modkit extract --reference $3 --threads 64 --queue-size 2000000 --ignore h -n 50000 --read-calls-path $2 --mapped-only --cpg $1 null
#Motif restrict to CpG sites, other motif could be specified as well
#read-calls to using the same thresholding algorithm as pile up
#mapped only, only extract the reads aligend to the reference genome
#run in 64 threads, 2000000 reads can be in memory at sametime