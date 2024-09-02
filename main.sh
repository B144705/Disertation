# This is not a script for run directly but a instruction about what to input and ouput should be for each step
#login to eddie
ssh UUN@eddie.ecdf.ed.ac.uk
#all the scripts run in main working directory
####################################################################################################
#EXTRACT THE METHYLATION INFORMATION INTO SINGLE READ LEVEL

#MAIN OUTPUTS
#bed file contianing the single read methylation information

#Define directories and variables
methylation_extract_script=# path to singleread_methylation_extract.sh
modified_bam=#path to modified bam file
tsv_out=#path to tsv file contian single read methylation information
ref_genome=#path to hg38 reference genome
bed_cov=#path to mod_tsv2bed.py
bed_single=#path to store bedfile contian single read methylation information
#TSV FILE contain methylation information
qsub ${methylation_extract_script} ${modified_bam} ${tsv_out} ${ref_genome}
#login to interactive node
#Convert tsv file to bed file contain single read information
qlogin -pe interactivemem 1 -l h_vmem=32G #require 32G RAM
module load python/3.11.4
python ${bed_cov} -i ${tsv_out} > ${bed_single}
####################################################################################################
#BULK METHYLATION OF INDIVIDUAL CpG SITE

#MAIN OUTPUTS
#bedmethyl file containing the methylation level for individual CpG site
#Define directories and variables
modified_bam=#path to modified bam file
ref_genome=#path to hg38 reference genome
bulk_meth=#path to store the bulk methylation bedfile
qlogin -pe interactivemem 32 -l h_vmem=2G # Request 16 core cpu
export PATH="MODKIT/PATH:$PATH" #load the modkit tool
modkit pileup --preset traditional --threads 32 --queue-size 1000000 --ref ${ref_genome} ${modified_bam} > ${bulk_out}
#The –preset traditional setting in modkit pileup was used so that:
# We only include reference CGs;
# Only C and 5mC are reported rather than other modifications;
# Data is aggregated across strands
####################################################################################################
#SINGLE-READ STATISTICS
# MAIN OUTPUTS:
# bed files containing read-level statistics
# bed files containing read-level statistics along with window info for all reads that align entirely within one window
#Define directories and variables
bed_single=#path to store bedfile contian single read methylation information
R_single_read_stats_script= # path to single_read_statistics.R R script
single_red_stat_dir=#path to directroy store single read stat
window_bed_file= # path to bed file containing genomic windows
dataset=#name of orginal dataset
qlogin -pe interactivemem 1 -l h_vmem=32G
module load R/4.3.0
Rscript ${R_single_read_stats_script} ${bed_single} ${single_red_stat_dir}/${dataset}
# Use bedtools to identify the reads that are in each genomic window
bedtools intersect -bed -wa -wb -a ${single_red_stat_dir}/${dataset} -b ${window_bed_file} -f 1 | awk 'OFS="\t" {print $1,$2,$3,$4,$8,$9,$10,$11,$12,$13,$14,$15,$19}' > singe_read_stat_with_window_info.bed
####################################################################################################
#GENOMIC-WINDOWS STATISTICS
# MAIN OUTPUTS
# tsv files containing the mean single read stat for each window

#Define directories and variables
window_stat_script=#path to calculate_mean_window_stats.R
single_readstat_windowinfo=#path to bed files containing read-level statistics along with window info for all reads that align entirely within one window
window_stat_dir=#path to directroy store genomic windows mean statistics
qlogin -pe interactivemem 1 -l h_vmem=32G
module load R/4.3.0
Rscript ${window_stat_script} ${single_readstat_windowinfo} ${window_stat_dir}/windowstat.bed

####################################################################################################
# MOST/LEAST HETEROGENEOUS WINDOWS (DEFINED VIA RTS)
# MAIN OUTPUTS
# bed files of 10% of windows with highest RTS and the 10% of windows with the lowest RTS (most heterogenous and least heterogenous)
#Define directories and variables
RTS_script=#path to highest_lowest_RTS_windows.sh
window_stat_tsv=#path to tsv file containing the mean single read stat for each window
window_stat_dir=#path to directroy store genomic windows mean statistics
mkdir ${window_stat_dir}/high_low_heterogenity_windows
qsub ${RTS_script} ${window_stat_tsv} ${window_stat_dir}/high_low_heterogenity_windows


####################################################################################################
# IDENTIFY PMDS AND NON-PMDS

# MAIN OUTPUTS
# bed files containing the genomic locations of PMDs and non-PMDs (each with min length 200kb)

# Define directories and variables

pmd_nonpmd_identification_script= # path to identify_PMDs_nonPMDs.sh script
main_dir= # path to main working directory
bulk_meth=#path store the bulk methylation bedfile
chrom_sizes= # path to sorted chromosome size file (sort -k1,1)

qsub ${pmd_nonpmd_identification_script} ${main_dir} ${CpG_bulk_meth_info} ${chrom_sizes}
####################################################################################################
#CACULATION OF FOLD ENRICHMENT IN PMD/non-PMD
# MAIN OUTPUTS
# Fold enrichment of most heterogenous and least heterogenous windows in PMD/non-PMD
#Define directories and variables
PMD_enrich_analysis=#path to PMD_enrich_analysis.sh
window_all_bed=# path to bed file containing ALL genomic windows
window_high_RTS_bed=#path to bed file containing windows with top 10% RTS
window_low_RTS_bed=#path to bed file contain windwos with bottom 10% RTS
PMD_info=#path to bed file contain PMD annotation
non_PMD_info=#path to bed file contain non-PMD anootation
enrich_dir=#directory store fold enrichment information

qlogin -pe interactivemem 1 -l h_vmem=32G
bash ${PMD_enrich_analysis} ${window_all_bed} ${window_high_RTS_bed} ${window_low_RTS_bed} ${PMD_info} ${non_PMD_info} ${enrich_dir}
####################################################################################################
#CACULATION OF FOLD ENRICHMENT IN ChromHMM states
# MAIN OUTPUTS
# FOLD enrichment of most heterogenous and least heterogenous windows in different ChromHMM states
#Define directories and variables
chromHMM_enrich=#path to chromHMM_enrichment.sh
window_all_bed=# path to bed file containing ALL genomic windows
window_high_RTS_bed=#path to bed file containing windows with top 10% RTS
window_low_RTS_bed=#path to bed file contain windwos with bottom 10% RTS
chromHMM_bed=#path to chromHMM bedfile
enrich_dir=#directory store fold enrichment information

qlogin -pe interactivemem 1 -l h_vmem=32G
bash ${chromHMM_enrich} ${window_all_bed} ${window_high_RTS_bed} ${window_low_RTS_bed} ${chromHMM_bed} ${enrich_dir}
####################################################################################################
#CACULATION OF SINGLE READ DISTANCE DEPENDENT CORRELATION
# MAIN OUPUT
# distancse dependet correlation of all the individual reads>=100CpG aligned to the region of intereset
#Step 1
#Extract the reads aligned to the regions of interest(>=100CpG)
#Define directories and variable
single_read_meth_stat=#path to single read methylation bedfile contain single read statistics
locations_bed=#path to bed file containing the genomic locations of interest

qlogin -pe interactivemem 1 -l h_vmem=32G
module load igmm/apps/BEDTools/2.30.0
bedtools intersect -f 1 -wa -a ${single_read_meth_stat} -b ${locations_bed} | awk '{if($8>=100) print$0}' > single_read_read_aligend_locations_100CpG.bed
#Step 2
#Extract the distance between CpG site within single read
#Define directories and variable
distance_extract_script=#path to methylation_call_extract.sh
read_aligend_interest=#path to single_read_read_aligend_locations_100CpG.bed
distance_info_dir=#path to directroy store distance information
combine_script=#path to distance_info_combine.sh
#OUTPUT OF STEP2
#read_distances.txt extract all the distance between neighborhood CpG site
#methcall_1.txt keep all the methylation call except last
#methcall_2.txt keep all the methylation call except first
#read_dist_info.bed combine above 3 to get the bed file contain distance inforamtion about CpG sites within each read
bash ${distance_extract_script} ${read_aligend_interest} ${distance_info_dir}
bash ${combine_script} ${read_aligend_interest} methcall_1.txt methcall_2.txt read_distances.txt
#Step 3
#Caculate the distance dependent correlation in batches
#Define directories and variable
distance_dependent_correlation_script=#path to autocorrealtion.sh
read_batch=#path to directiroy contain subdirectory of each batches
batch_txt=#path to batch.txt contain the list of the name for all the batches
correlation_dir=#path to directory store the caculated distance dependent correlation
n=#number of batches
#spilit the large dist_info file into small batch each contain 10000 reads
split -l 10000 --additional-suffix=.bed -d ${100CpG_read_dist_info.bed} ${read_batch}/100CpG_read_dist_info
qsub -t 1-$n ${distance_dependent_correlation_script} ${batch_txt} ${read_batch} ${correlation_dir}







