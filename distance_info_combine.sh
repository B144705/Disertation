#$1 single_read_meth_100CpG.bed
#$2 100CpG_methcall_1.txt
#$3 100CpG_methcall_2.txt
#$4 100CpG_read_distances.txt

awk 'OFS="\t" {print $1, $2, $3, $4}' $1 | paste -d '\t' - $2 $3 $4 | awk -F '\t' '{num_distances=split($5,a,","); num_call=split($7,a,","); if(num_distances==num_call) {print $0}}' > 100CpG_read_dist_info.bed
