PROMOTEROME_PARENT_DIR="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome"
REP="rep2" # change
PROMOTEROME_DIR="Jurkat_18hr_"$REP # change
FASTQ_R1="Jurkat_18hr_"$REP"_R1.fastq.gz" # change
FASTQ_R2="Jurkat_18hr_"$REP"_R2.fastq.gz" # change
RAWDATA_DIR="/mnt/isilon/sfgi/rawData/wells/captureC/Jurkat/Jurkat_18hr/Jurkat_18hr_"$REP  # change

############## CREATE ANALYSIS DIR FOR EACH REP ################
mkdir $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/
mkdir $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

### CREATE FASTQ SYMBOLIC LINK FOR EACH REP
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR
ln -s $RAWDATA_DIR/$FASTQ_R1
ln -s $RAWDATA_DIR/$FASTQ_R2

### PREPARE FOR HICUP
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup

#### CREATE FASTQ CHUNK FOR R1 AND R2
sbatch -o split_1_1.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R1 200000000 1 R1
sbatch -o split_1_2.log ~/captureC/scripts/bash/splitFastqGz.sh ../$FASTQ_R2 200000000 1 R2

### AFTER CREAT CHUNK CHECK WHETHER TOTAL IS RIGHT
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
total_split=$(ls lane_1_R1* | wc -l)
total_split_0base=$((total_split-1))
files=($(ls -lhtr lane_1_R1* | awk '{print $9}'))
last_i=$((${#files[@]} - 1))
reminder=$(wc -l ${files[$last_i]} | awk '{num=$1/4; print num}')
total=$(cat $RAWDATA_DIR/readCount.txt)
cal_total=$((total_split_0base*50000000+reminder))
if [ $total = $cal_total ]; then
        echo "$cal_total is correct"
else
        echo "splitted read number not match"
fi

total_split=$(ls lane_1_R2* | wc -l)
total_split_0base=$((total_split-1))
files=($(ls -lhtr lane_1_R2* | awk '{print $9}'))
last_i=$((${#files[@]} - 1))
reminder=$(wc -l ${files[$last_i]} | awk '{num=$1/4; print num}')
total=$(cat $RAWDATA_DIR/readCount.txt)
cal_total=$((total_split_0base*50000000+reminder))
if [ $total = $cal_total ]; then
        echo "$cal_total is correct"
else
        echo "splitted read number not match"
fi

### CHANGE FILE NAME
files=($(ls -lhtr lane_1_R1* | awk '{print $9}'))
for i in `seq 1 ${#files[@]}`;do 
        j=$((i-1))
        newfile="lane_1_"$i"_R1.fq"
        if [[ $newfile != $files[$j] ]]; then
                # echo "$newfile ${files[$j]}"
                mv ${files[$j]} $newfile
        fi
done

files=($(ls -lhtr lane_1_R2* | awk '{print $9}'))
for i in `seq 1 ${#files[@]}`;do 
        j=$((i-1))
        newfile="lane_1_"$i"_R2.fq"
        if [[ $newfile != $files[$j] ]]; then
                # echo "$newfile ${files[$j]}"
                mv ${files[$j]} $newfile
        fi
done

### FOR EACH LANE, CREATE CHUNK DIRECTORIES
for i in `seq 1 $total_split`;do mkdir lane_1_${i};done

### GZIP FASTQ
for i in `seq 1 $total_split`;do sbatch -o logs -J "gzip_1_"${i}"_R1"  ~/captureC/scripts/bash/gzip_mv.sh "lane_1_"$i"_R1.fq";done
for i in `seq 1 $total_split`;do sbatch -o logs -J "gzip_1_"${i}"_R2"  ~/captureC/scripts/bash/gzip_mv.sh "lane_1_"$i"_R2.fq";done


### PREPARE JOB BASH AND CONFIG FILE FOR HICUP (hg19) - can run hg19|hg38 and Arima|DpnII
# more info `perl ~/captureC/scripts/perl/writeHicupScripts.v2.pl`
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
total_split=$(ls -d lane_1_*/ | wc -l)
perl ~/captureC/scripts/perl/writeHicupScripts.v2.pl --genome hg19 --enzyme DpnII $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup $total_split 1 8G

for i in `seq 1 $total_split`;do 
        cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/lane_1_${i}
        sbatch -t 48:00:00 -c 2 --mem-per-cpu 8G -o hicup.log -J "hicup_"$i hicup.sh
        cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
done



### CHECK HICUP RUNNING (IN DIFFERENT CHUNK) ARE SUCCESSFULLY COMPLETED
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
grep "HiCUP processing complete" lane*/hicup.log | wc -l
# check whether it matches split chunk number
total_split=$(ls -d lane_1_*/ | wc -l)
echo $total_split


echo "source ~/.bashrc
module load R/current
files=(\$(ls $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/*/*.hicup.bam))
perl ~/captureC/scripts/perl/hicupMergeDedup.v0.7.4.pl --threads 4 \${files[@]}
" > hicupMergeDedup.sh
sed -i '1i#!/bin/bash' hicupMergeDedup.sh
# qsub -cwd -l h_vmem=4G -pe smp 4 -N hicupMergeDedup -o hicupMergeDedup.log hicupMergeDedup.sh
sbatch -t 48:00:00 --mem-per-cpu 4G -c 4 -J hicupMergeDedup hicupMergeDedup.sh

### CAPUTRE
# CONVERT TO CHICAGO INPUT (1 FRAGMENT RESOLUTION)
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
DIR_1FRAG=$REP"_1frag"

bash ~/captureC/scripts/bash/run_bam2chicago.sh \
-k DpnII -g hg19 -f 1frag \
-m 64G -s 1 \
$DIR_1FRAG $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam

# CONVERT TO CHICAGO INPUT (4 FRAGMENT RESOLUTION)
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
DIR_1FRAG=$REP"_4frag"

bash ~/captureC/scripts/bash/run_bam2chicago.sh \
-k DpnII -g hg19 -f 4frag \
-m 64G -s 1 \
$DIR_1FRAG $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/merged.dedup.bam


### COLLECT INFO AND OUTPUT AS HICUP SUMMARY
cd $PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup
total_split=$(ls -d lane_1_*/ | wc -l)
deduplicator_FILE=$(ls hicup_deduplicator_summary*)
SUMMARY_FILE=$PROMOTEROME_DIR"_hicupSummary.txt"

DIR_1FRAG=$REP"_1frag" # optional only used for capture C
chinput_FILE=$DIR_1FRAG".chinput" # optional only used for capture C


perl ~/captureC/scripts/perl/collectHicupResults.v2.pl 1 $total_split \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$deduplicator_FILE \
$PROMOTEROME_PARENT_DIR/$PROMOTEROME_DIR/hicup/$DIR_1FRAG/$chinput_FILE \
> $SUMMARY_FILE


DIR_4FRAG=$REP"_4frag" # optional only used for capture C
chinput_FILE=$DIR_4FRAG".chinput" # optional only used for capture C
TAG_NUM=$(sed '6q;d' $SUMMARY_FILE | awk '{print $24}')
perl ~/captureC/scripts/perl/captureStatsFromChinput.pl $DIR_4FRAG/$chinput_FILE $TAG_NUM > 4frag_captured.txt

bash ~/captureC/scripts/bash/summarize_hicup.v5.sh -o $PROMOTEROME_DIR"_ShortSummary.txt" $PROMOTEROME_DIR"_hicupSummary.txt"

# CLEAN-UP (kept filt.bam hicup.bam)
rm -f lane_*_*/*fastq.gz
rm -f lane_*/*fq.gz
rm -f lane_*/*pair.bam
rm merged.bam