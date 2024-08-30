#!/bin/bash

# Specifications for eddie

#$ -N window_mean_stats_job                     # Name the job
#$ -pe sharedmem 1                              # Request 1 core cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=32G                                 # Request 32G memory per cpu
#$ -cwd
#$ -l h_rt=48:00:00                              # Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh
# Arguments:
# $1: path to single-molecule window stastics bedfile
# $2: path to bed file with genomic windows


# Load the required modules

module load igmm/apps/BEDTools/2.30.0

# Extract locations of all windows for which there was sufficient data to perform the single-molecule analysis
# Also save the number of these windows as a variable
name=$(basename $1 _100KB)
awk '{print $1}' $1 | grep -wf - $2 > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_analysed.bed
num_windows=$(wc -l nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_analysed.bed | awk '{print $1}')
# Create bedgraphs with the mean methylation and RTS for each window
sort -k1,1 $1 | join -1 1 -2 4 - <(sort -k4,4 nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_analysed.bed) | awk 'OFS="\t" {print $8,$9,$10,$2}' > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_mean_meth.bedGraph
sort -k1,1 $1 | join -1 1 -2 4 - <(sort -k4,4 nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_analysed.bed) | awk 'OFS="\t" {print $8,$9,$10,$6}' > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_windows_RTS.bedGraph


# Extract the 10% of windows with highest RTS and the 10% of windows with the lowest RTS
# Extract locations of these windows and save as bed files

awk 'NR==1; NR > 1 {print $0 | "sort -nr -k6,6"}' $1 | head -$(((${num_windows}+10-1)/10)) > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_top_10_percent_RTS.tsv
awk 'NR==1; NR > 1 {print $0 | "sort -n -k6,6"}' $1 | head -$(((${num_windows}+10-1)/10)) > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_bottom_10_percent_RTS.tsv

awk '{print $1}' nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_top_10_percent_RTS.tsv | grep -wf - $2 > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_top_10_percent_RTS_windows.bed
awk '{print $1}' nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_bottom_10_percent_RTS.tsv | grep -wf - $2 > nanopore_dnmt1i_2024/100KB_analysis/high_low_heterogenity_windows/${name}_bottom_10_percent_RTS_windows.bed






