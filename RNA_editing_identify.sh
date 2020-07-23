#!/usr/bin/bash
for i in {1..2}
do
###  #First mapping ##########
bwa aln -t 20 "$1" "$2"_"$i".clean.fq.gz -f "$2"_"$i".sai
bwa samse -n 50 "$1" "$2"_"$i".sai "$2"_"$i".clean.fq.gz  -f "$2"_"$i".aln.sam
samtools view -b -f4 "$2"_"$i".aln.sam -o "$2"_"$i".aln.um.bam
bam2fastx -o "$2"_"$i".aln.um.fastq "$2"_"$i".aln.um.bam
bwa mem -M -t 20 -k 50 "$1" "$2"_"$i".aln.um.fastq > "$2"_"$i".mem.sam 
samtools view -b -f4 "$2"_"$i".mem.sam -o "$2"_"$i".mem.um.bam
bam2fastx -o "$2"_"$i".mem.um.fastq "$2"_"$i".mem.um.bam
rm "$2"_"$i".sai  "$2"_"$i".mem.um.bam

perl Fastq_transform.pl "$2"_"$i".mem.um.fastq

#a2g++
bwa aln -t 5 -n 2 -o 0 -N "$1".a2g a2g."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2g temp.sai a2g."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2g++.origReads.sam
perl filter_sam.pl a2g++.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2g++.filterReads.sam
rm tem* a2g++.origReads.sam
#a2g+-
bwa aln -t 5 -n 2 -o 0 -N "$1".a2g t2c."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2g temp.sai t2c."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2g+-.origReads.sam
perl filter_sam.pl a2g+-.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2g+-.filterReads.sam
rm tem* a2g++.origReads.sam
#a2g-+
bwa aln -t 5 -n 2 -o 0 -N "$1".t2c a2g."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".t2c temp.sai a2g."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2g-+.origReads.sam
perl filter_sam.pl a2g-+.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2g-+.filterReads.sam
rm tem* a2g++.origReads.sam
#a2g--
bwa aln -t 5 -n 2 -o 0 -N "$1".t2c t2c."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".t2c temp.sai t2c."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2g--.origReads.sam
perl filter_sam.pl a2g--.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2g--.filterReads.sam
rm tem* a2g++.origReads.sam
#a2c++
bwa aln -t 5 -n 2 -o 0 -N "$1".a2c a2c."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2c temp.sai a2c."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2c++.origReads.sam
perl filter_sam.pl a2c++.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2c++.filterReads.sam
rm tem* a2g++.origReads.sam
#a2c+-
bwa aln -t 5 -n 2 -o 0 -N "$1".a2c t2g."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2c temp.sai t2g."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2c+-.origReads.sam
perl filter_sam.pl a2c+-.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2c+-.filterReads.sam
rm tem* a2g++.origReads.sam
#a2c-+
bwa aln -t 5 -n 2 -o 0 -N "$1".t2g a2c."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".t2g temp.sai a2c."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2c-+.origReads.sam
perl filter_sam.pl a2c-+.origReads.sam unique_simple_repeats.txt >"$2"_"$i".a2c-+.filterReads.sam
rm tem* a2g++.origReads.sam
#a2c--
bwa aln -t 5 -n 2 -o 0 -N "$1".t2g t2g."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".t2g temp.sai t2g."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2c--.origReads.sam
perl filter_sam.pl a2c--.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2c--.filterReads.sam
rm tem* a2g++.origReads.sam
#g2c++
bwa aln -t 5 -n 2 -o 0 -N "$1".g2c g2c."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".g2c temp.sai g2c."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > g2c++.origReads.sam
perl filter_sam.pl g2c++.origReads.sam unique_simple_repeats.txt > "$2"_"$i".g2c++.filterReads.sam
rm tem* a2g++.origReads.sam
#g2c+-
bwa aln -t 5 -n 2 -o 0 -N "$1".g2c c2g."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".g2c temp.sai c2g."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > g2c+-.origReads.sam
perl filter_sam.pl g2c+-.origReads.sam unique_simple_repeats.txt > "$2"_"$i".g2c+-.filterReads.sam
rm tem* a2g++.origReads.sam
#a2t++
bwa aln -t 5 -n 2 -o 0 -N "$1".a2t a2t."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2t temp.sai a2t."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2t++.origReads.sam
perl filter_sam.pl a2t++.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2t++.filterReads.sam
rm tem* a2g++.origReads.sam
#a2t+-
bwa aln -t 5 -n 2 -o 0 -N "$1".a2t t2a."$2"_"$i".mem.um.fastq -f temp.sai 
bwa samse -n 50 "$1".a2t temp.sai t2a."$2"_"$i".mem.um.fastq > temp.sam 
samtools view -F 4 -S -h -o temp3.sam temp.sam 
perl Get_orig_read.pl temp3.sam "$2"_"$i".clean.fq.gz > a2t+-.origReads.sam
perl filter_sam.pl a2t+-.origReads.sam unique_simple_repeats.txt > "$2"_"$i".a2t+-.filterReads.sam
rm tem* a2g++.origReads.sam

for bamFILE in "$2"_"$i".*filterReads.sam;
do
perl analyse_mm.pl "$1".fa $bamFILE >> "$2"_"$i".mem.um.temp_file
done
perl sort_R_read.pl "$2"_"$i".mem.um.temp_file "$2"_"$i".mem.um.temp_statistic > "$2"_"$i".mem.um.analyseMM
rm "$2"_"$i".mem.um.temp_file 
done
####Find RNA editing######
perl Detect_editing.pl "$2"_1.mem.um.analyseMM "$2"_1_results "$2"_2.aln.sam "$2"_2.mem.sam
perl Detect_editing.pl "$2"_2.mem.um.analyseMM "$2"_2_results "$2"_1.aln.sam "$2"_1.mem.sam
rm "$2"_1.aln.sam "$2"_1.mem.sam "$2"_2.aln.sam "$2"_2.mem.sam

