### CHROMHMM ENRICHMENTS ###
# $1: bed file containing all windows
# $2: bed file containing lowest RTS windows
# $3: bed file containig highest RTS windows
# $4: chromHMM bed file
# $5: path to output folder

# Calculate the number of bases in all windows
all_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $1)

# Calculate the number of bases in the lowest correlation windows
lowest_RTS_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $2)

# Calculate the number of bases in the highest correlation windows
highest_RTS_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3)


# For each chromatin state extract the overlap and proportion overlap with all windows, low RTS windows and high RTS windows
module load igmm/apps/BEDTools/2.30.0
mkdir $5/specific_state_overlaps
for i in $(awk '{print $4}' $4 | sort -k1,1 | uniq | cat)
do
 grep -w $i $4 | bedtools intersect -wao -a $1 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $5/specific_state_overlaps/window_proportion_overlap_$i.txt.gz;
 grep -w $i $4 | bedtools intersect -wao -a $2 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $5/specific_state_overlaps/bottom_10_percent_RTS_proportion_overlap_$i.txt.gz;
 grep -w $i $4 | bedtools intersect -wao -a $3 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $5/specific_state_overlaps/top_10_percent_RTS_proportion_overlap_$i.txt.gz;
done

# Extract the overlap between all windows and the chromHMM annotations in one file
bedtools intersect -a $4 -b $1 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $5/window_overlap_annotations.bed

# For all the windows analysed, calculate the % of bp cooresponding to each annotation

rm $5/all_window_overlap_chromHMM_percents.txt;
touch $5/all_window_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $4 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $5/window_overlap_annotations.bed | awk -v a="${all_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $5/all_window_overlap_chromHMM_percents.txt;
 done



# Repeat for the lowest RTS windows

# Extract the overlap between the lowest RTS windows and the chromHMM annotations
bedtools intersect -a $4 -b $2 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $5/bottom_10_percent_RTS_overlap_annotations.bed

# For the lowest RTS windows, calculate the % of bp cooresponding to each annotation
rm $5/bottom_10_percent_RTS_overlap_chromHMM_percents.txt
touch $5/bottom_10_percent_RTS_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $4 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $5/bottom_10_percent_RTS_overlap_annotations.bed | awk -v a="${lowest_RTS_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $5/bottom_10_percent_RTS_overlap_chromHMM_percents.txt;
 done



# Repeat for the highest RTS windows

# Extract the overlap between the highest RTS windows and the chromHMM annotations
bedtools intersect -a $4 -b $3 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $5/top_10_percent_RTS_overlap_annotations.bed

# For the highest RTSs windows, calculate the % of bp cooresponding to each annotation
rm $5/top_10_percent_RTS_overlap_chromHMM_percents.txt
touch $5/top_10_percent_RTS_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $4 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $5/top_10_percent_RTS_overlap_annotations.bed | awk -v a="${highest_RTS_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $5/top_10_percent_RTS_overlap_chromHMM_percents.txt;
 done



# NOTE THAT ENRICHMENTS CAN THEN BE OBTAINED BY DIVIDING VALUES IN bottom_10_percent_RTS_overlap_chromHMM_percents.txt
# and top_10_percent_RTS_overlap_chromHMM_percents.txt BY VALUES IN all_windows_percent_RTS_overlap_chromHMM_percents.txt.
# THIS IS NOT DONE HERE AS INSTEAD MATHEMATICA WAS USED TO OBTAIN THE ENRICHMENTS.