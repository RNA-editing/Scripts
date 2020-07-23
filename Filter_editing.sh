#!/bin/bash
# echo "Welcome using RNA editing filter!"

echo "Fiter Reads for $1" 
awk '{if($NF<=($(NF-1)+$NF)*0.2 && $(NF-1)>=3)print}' "$1"-1.UE.list> "$1"-1.UE.filter.list
awk '{if($NF<=($(NF-1)+$NF)*0.2 && $(NF-1)>=3)print}' "$1"-2.UE.list> "$1"-2.UE.filter.list
#wc -l "$i"_1.UE.filter.list>> Hyper_editing.reads.number.txt
#wc -l "$i"_2.UE.filter.list>> Hyper_editing.reads.number.txt
grep 'Aligns to:' "$1"-1.UE.Details|awk  '{print $NF"\t"$(NF-3)"\t"$(NF-6)}'|sed 's/,//g'|awk -F '[()]' '{print $0"\t"$2}'|awk '{if($4=="-")$4="0";else if($4="+")$4="1"; print}' > "$1"-1.position.txt
grep 'Edit indexes:' "$1"-1.UE.Details|awk '{print $3}'|sed 's/,//g'> "$1"-1.ps.txt
paste "$1"-1.position.txt "$1"-1.ps.txt> "$1"-1.pos.txt
perl extract_Reads.pl "$1"-1.UE.filter.list "$1"-1.pos.txt >"$1"-1.cluster.txt
perl SNPs_splicing_filter.pl "$1"-1.cluster.txt "$2" "$3"> "$1"-1.site.txt
rm "$1"-1.pos.txt "$1"-1.ps.txt "$1"-1.position.txt
grep 'Aligns to:' "$1"-2.UE.Details|awk  '{print $NF"\t"$(NF-3)"\t"$(NF-6)}'|sed 's/,//g'|awk -F '[()]' '{print $0"\t"$2}'|awk '{if($4=="-")$4="0";else if($4="+")$4="1"; print}' >"$1"-2.position.txt
grep 'Edit indexes:' "$1"-2.UE.Details|awk '{print $3}'|sed 's/,//g'> "$1"-2.ps.txt
paste "$1"-2.position.txt "$1"-2.ps.txt> "$1"-2.pos.txt
perl Cluster.pl "$1"-2.UE.filter.list "$1"-2.pos.txt >"$1"-2.cluster.txt
perl SNPs_splicing_filter.pl "$1"-2.cluster.txt "$2" "$3"> "$1"-2.site.txt
rm "$1"-2.pos.txt "$1"-2.ps.txt "$1"-2.position.txt
###merge two reads
echo "Merge Pair Reads" 
awk '{if($4==0&&$3=="A2G")$3="T2C";else if($4==0&&$3=="A2T")$3="T2A";else if($4==0&&$3=="A2C")$3="T2G";else if($4==0&&$3=="T2G")$3="A2C";else if($4==0&&$3=="T2C")$3="A2G"; else if($4==0&&$3=="T2A")$3="A2T";else if($4==0&&$3=="C2A")$3="G2A";else if($4==0&&$3=="C2G")$3="T2C";else if($4==0&&$3=="C2T")$3="G2A";else if($4==0&&$3=="G2A")$3="C2T";else if($4==0&&$3=="G2T")$3="C2A";else if($4==0&&$3=="G2C")$3="C2G"; print}' "$1"-2.site.txt|sort -k1,1 -k2,2n> "$1"-2.position
awk '{if($4==1&&$3=="A2G")$3="T2C";else if($4==1&&$3=="A2T")$3="T2A";else if($4==1&&$3=="A2C")$3="T2G";else if($4==1&&$3=="T2G")$3="A2C";else if($4==1&&$3=="T2C")$3="A2G"; else if($4==0&&$3=="T2A")$3="A2T";else if($4==1&&$3=="C2A")$3="G2A";else if($4==1&&$3=="C2G")$3="T2C";else if($4==1&&$3=="C2T")$3="G2A";else if($4==1&&$3=="G2A")$3="C2T";else if($4==1&&$3=="G2T")$3="C2A";else if($4==1&&$3=="G2C")$3="C2G"; print}' "$1"-1.site.txt|sort -k1,1 -k2,2n|uniq> "$1"-1.position
awk '{if($3==0&&$2=="A2G")$2="T2C";else if($3==0&&$2=="A2T")$2="T2A";else if($3==0&&$2=="A2C")$2="T2G";else if($3==0&&$2=="T2G")$2="A2C";else if($3==0&&$2=="T2C")$2="A2G"; else if($3==0&&$2=="T2A")$2="A2T";else if($3==0&&$2=="C2A")$2="G2A";else if($3==0&&$2=="C2G")$2="T2C";else if($3==0&&$2=="C2T")$2="G2A";else if($3==0&&$2=="G2A")$2="C2T";else if($3==0&&$2=="G2T")$2="C2A";else if($3==0&&$2=="G2C")$2="C2G"; print}' "$1"-2.cluster.txt> "$1"-2.reads.position
awk '{if($3==1&&$2=="A2G")$2="T2C";else if($3==1&&$2=="A2T")$2="T2A";else if($3==1&&$2=="A2C")$2="T2G";else if($3==1&&$2=="T2G")$2="A2C";else if($3==1&&$2=="T2C")$2="A2G"; else if($3==0&&$2=="T2A")$2="A2T";else if($3==1&&$2=="C2A")$2="G2A";else if($3==1&&$2=="C2G")$2="T2C";else if($3==1&&$2=="C2T")$2="G2A";else if($3==1&&$2=="G2A")$2="C2T";else if($3==1&&$2=="G2T")$2="C2A";else if($3==1&&$2=="G2C")$2="C2G"; print}' "$1"-1.cluster.txt> "$1"-1.reads.position
echo -e -n "Chr\tPosition\tEditing_type\tStrand\tReads_postion\n" >"$1".RNA_editing_events.txt
cat "$1"-2.position "$1"-1.position|sort -k1,1 -k2,2n>>"$1".RNA_editing_events.txt
echo -e -n "Chr\tPosition\tEditing_type\n" >"$1".uniq.sites.txt
cat "$1"-2.position "$1"-1.position|awk '{print $1"\t"$2"\t"$3}'|uniq>>"$1".uniq.sites.txt
echo -e -n "Reads_position\tEditing_type\tStrand\tEditing_Position_on_reads" >"$1".hyper_reads.txt
cat "$1"-1.reads.position "$1"-2.reads.position|sort>>"$1".hyper_reads.txt
#awk '{print $1"_"$2"_"$3}' "$1".uniq.position >"$1".uniq.id
rm "$1"-1.reads.position
rm "$1"-2.reads.position
rm "$1"-1.cluster.txt "$1"-2.cluster.txt "$1"-1.site.txt "$1"-2.site.txt "$1"-1.position "$1"-2.position
#rm "$1".uniq.position
for j in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
do 
perl Cluster.pl "$1".uniq.sites.txt "$j">> "$1".cluster.tem.position
done
echo -e -n "Chr\tEditing_type\tEditing_Position\n" >"$1".cluster.txt
grep -v '^$' "$1".cluster.tem.position>>"$1".cluster.txt
rm "$1".cluster.tem.position "$1"-1.UE.filter.list  "$1"-2.UE.filter.list 
#awk -F '[\t,]' '{print $1"\t"$3"\t"$(NF-1)"\t"$2}' "$1".cluster.position > "$1".cluster.bed

#####
#for j in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
#do 
#echo "Merge cluster for $j"
##for i in $x
##do
#grep ''"$j"'' "$1".cluster.bed|awk '{print $1"\t"$2"\t"$3"\t"'"$i"'}' > "$1"."$j".bed
##done
#cat *."$j".bed|sort -k1,1 -k2,2n > All."$j".bed
#bedtools merge -i All."$j".bed -c 4 -o collapse >All.merge."$j".bed
#rm *[0-9]."$j".bed
#rm All."$j".bed
#done
##### Count Cluster type
#echo "Count type"
#echo -e -n "Type" >>Cluster_type.txt
#echo -e -n "Type" >>Uniq_sites_type.txt
#echo -e -n "Type" >>Editing_reads_type.txt
#echo -e -n "Type" >>Event_reads_type.txt
##for j in $x
##do
#echo -e -n "\t$j"  >>Cluster_type.txt
#echo -e -n "\t$j"  >>Uniq_sites_type.txt
#echo -e -n "\t$j"  >>Editing_reads_type.txt
#echo -e -n "\t$j"  >>Event_reads_type.txt
##done
#echo ""  >>Cluster_type.txt
#echo ""  >>Uniq_sites_type.txt
#echo ""  >>Editing_reads_type.txt
#echo ""  >>Event_reads_type.txt
#for i in A2G A2C A2T C2G C2A C2T G2C G2A G2T T2A T2G T2C
#do
#echo -n "$i"  >>Cluster_type.txt
#echo -n "$i"  >>Uniq_sites_type.txt
#echo -n "$i"  >>Editing_reads_type.txt
#echo -n "$i"  >>Event_reads_type.txt
##for j in $x 
##do
#echo -e -n "\t"  >>Cluster_type.txt
#echo -e -n "\t"  >>Uniq_sites_type.txt
#echo -e -n "\t"  >>Editing_reads_type.txt
#echo -e -n "\t" >>Event_reads_type.txt
#grep ''"$i"'' SRR"$1".cluster.position|wc -l|tr -d '\n'  >>Cluster_type.txt
#grep ''"$i"'' SRR"$1".uniq.position|wc -l|tr -d '\n'  >>Uniq_sites_type.txt
#grep ''"$i"'' SRR"$1".reads.position|wc -l|tr -d '\n'  >>Editing_reads_type.txt
#grep ''"$i"'' SRR"$1".event.site|wc -l|tr -d '\n'  >>Event_reads_type.txt
#done 
#echo "" >>Cluster_type.txt
#echo "" >>Uniq_sites_type.txt
#echo "" >>Editing_reads_type.txt
#echo "" >>Event_reads_type.txt
##done
#
#echo "All Finish...Enjoy!"
#
