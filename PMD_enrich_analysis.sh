### PMD/non-PMD ENRICHMENTS ###
#Arguments
#$1 bed file contain all windoes (autosome)
#$2 bed file contian windows with top 10% RTS
#$3 bed file contain windwos with bottom 10% RTS
#$4 bed file contian PMD informaion of corresponding cell line
#$5 bed file contian nonPMD informaion of corresponding cell line
#$6 Directory store PMD enrichment information 

module load igmm/apps/BEDTools/2.30.0
# Calculate the number of bases in all windows
all_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $1)


# Calculate the number of bases in the highest RTS windows
highest_RTS_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $2)

# Calculate the number of bases in the lowest RTS windows
lowest_RTS_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3)

# Extract the overlap and proportion overlap of all windows with PMDs and non-PMDs

all_pmd_overlap=$(bedtools intersect -a $4 -b $1 | awk '{Tot+=$3-$2} END {print Tot}')
all_nonpmd_overlap=$(bedtools intersect -a $5 -b $1 | awk '{Tot+=$3-$2} END {print Tot}')
all_pmd_prop_overlap=$(echo "${all_pmd_overlap}/${all_bases}" | bc -l)
all_nonpmd_prop_overlap=$(echo "${all_nonpmd_overlap}/${all_bases}" | bc -l)

# Extract the overlap and proportion overlap of the lowest RTS windows with PMDs and non-PMDs

lowest_RTS_pmd_overlap=$(bedtools intersect -a $4 -b $3 | awk '{Tot+=$3-$2} END {print Tot}') 
lowest_RTS_nonpmd_overlap=$(bedtools intersect -a $5 -b $3 | awk '{Tot+=$3-$2} END {print Tot}')

lowest_RTS_pmd_prop_overlap=$(echo "${lowest_RTS_pmd_overlap}/${lowest_RTS_bases}" | bc -l)
lowest_RTS_nonpmd_prop_overlap=$(echo "${lowest_RTS_nonpmd_overlap}/${lowest_RTS_bases}" | bc -l)

# Extract the overlap and proportion overlap of the highest RTS windows with PMDs and non-PMDs

highest_RTS_pmd_overlap=$(bedtools intersect -a $4 -b $2 | awk '{Tot+=$3-$2} END {print Tot}') 
highest_RTS_nonpmd_overlap=$(bedtools intersect -a $5 -b $2 | awk '{Tot+=$3-$2} END {print Tot}')
highest_RTS_pmd_prop_overlap=$(echo "${highest_RTS_pmd_overlap}/${highest_RTS_bases}" | bc -l)
highest_RTS_nonpmd_prop_overlap=$(echo "${highest_RTS_nonpmd_overlap}/${highest_RTS_bases}" | bc -l)

# Calculate enrichments and save to file

echo "PMD enrichment: $(echo "${lowest_RTS_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $6/low_RTS_pmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${lowest_RTS_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $6/low_RTS_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${highest_RTS_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $6/high_RTS_pmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${highest_RTS_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $6/high_RTS_nonpmd_enrichments.txt