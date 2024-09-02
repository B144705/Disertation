#!/bin/bash

# Specifications for eddie
#$ -N correlation_average			# Name the job
#$ -R y                             # Reserve requested nodes as they become available
#$ -l h_vmem=32G					# Request 32G memory
#$ -cwd
#$ -l h_rt=24:00:00					# Allocate a maximum of 24 hours to the job

#arguments
PMD_corr_output_dir="/gpfs/igmmfs01/eddie/sproul-lab/Gawain/3BKO_subPMD/D3BKO_K27_pmd_autocorrelation_reads" # path to output directory for correlations associated with PMDs reads
# Load bedtools
module load igmm/apps/BEDTools/2.30.0  
PMD_corr_dir="/gpfs/igmmfs01/eddie/sproul-lab/Gawain/3BKO_subPMD/D3BKO_K27_pmd_mean" # path to directory to store overall high correlation read outputs
for i in {0..28};\
 do echo $i; awk 'FNR>1' ${PMD_corr_output_dir}/D3BKO_K27_pmd_100CpG_WT_PMD_read_dist_info$i.bed/*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2,3 -o mean > ${PMD_corr_dir}/D3BKO_K27_pmd_100CpG_WT_PMD_read_dist_info_$i.tsv;
done
cat ${PMD_corr_dir}//D3BKO_K27_pmd_100CpG_WT_PMD_read_dist_info_*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2,3 -o mean > ${PMD_corr_dir}/all_vs_all_D3BKO_K27_pmd_100CpG_mean_distance_dependent_correlations.tsv
