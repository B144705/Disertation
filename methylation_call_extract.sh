#this script is use to extract the methylation status from single readbed fiels
#INPUT example: 0m2u3m4u5m17m > methcall_1: m,u,m,u,m methcall_2: u,m,u,m,m
#$1 input file
#extrac the distance information
awk '{print $5}' $1 | sed 's/m/,/g' - | sed 's/u/,/g' | sed 's/x/,/g'| sed 's/.$//' | sed 's/^.//' | sed 's/^.//' | sed 's/^$/NA/g' > read_distances.txt
#replace any nummber with comma,remove the first comma, replace empty string with NA
awk '{print $5}' $1 | sed 's/[0-9]\+/,/g' | sed 's/^.//' | sed 's/^$/NA/g' > temp.txt
#keep all the methylation call except lasst
awk -F ',' 'OFS="," {if (NF>1) {$NF=""; print $0} else {print "NA"}}' temp.txt | sed 's/,$//' > methcall_1.txt
#keep all the methylation call except first
awk -F ',' 'OFS="," {if (NF>1) {$1=""; print $0} else {print "NA"}}' temp.txt | sed 's/^,//' > methcall_2.txt
rm temp.txt