atac_dir="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT"
cd $atac_dir

mkdir -p $atac_dir/peaks

mkdir -p $atac_dir/peaks/individualRep
cd $atac_dir/peaks
files=($(ls ../*/out/peak/macs2/rep*/*.300K.filt.narrowPeak.gz))
for file in ${files[@]}; do
        prefix=$(echo $file | cut -d "/" -f 7 | cut -d "." -f 1 | cut -d "_" -f 1-3)
        newfile=$prefix".bed"
        zcat $file | cut -f 1,2,3 | sort -u | sort-bed - | awk -v prefix=$prefix -v i=1 'BEGIN {OFS="\t"}{name=prefix"_"i; i++; print $1,$2,$3,name}' > individualRep/$newfile
        peakNum=$(wc -l individualRep/$newfile | awk '{print $1}')
        echo -e "$newfile\t$peakNum"
done > individualRep/IndividualRepPeakNum.txt

####################### pooledRep 
mkdir -p $atac_dir/peaks/pooledRep
cd $atac_dir/peaks

files=($(ls ../*/out/peak/macs2/pooled_rep/*.300K.filt.narrowPeak.gz))
for file in ${files[@]}; do
        newfile=$(echo $file | cut -d "/" -f 7 | cut -d "." -f 1 | cut -d "_" -f 1-2)
        prefix=$(echo $newfile | sed "s/_//g")
        newfile=$newfile"_pooledRep.bed"
        zcat $file | cut -f 1,2,3 | sort -u | sort-bed - | awk -v p=$prefix 'BEGIN{OFS="\t"}{name=p"_"NR; print $0,name}'> pooledRep/$newfile
        peakNum=$(wc -l pooledRep/$newfile | awk '{print $1}')
        echo -e "$newfile\t$peakNum"
done > pooledRep/pooledRepPeakNum.txt

# filter pooledRep peaks by ind rep peaks
mkdir -p $atac_dir/peaks/pooledRep/compInd
cd $atac_dir/peaks
pooled_files=($(ls pooledRep/*_pooledRep.bed))
for pooled_file in ${pooled_files[@]}; do
        prefix=$(basename $pooled_file | cut -d "_" -f 1-2)
        ind_files=($(ls individualRep/$prefix*))
        for ind_file in ${ind_files[@]}; do
                ind_prefix=$(basename $ind_file | cut -d "." -f 1)
                # echo -e "$pooled_file\t$ind_file\t$ind_prefix"
                bedtools intersect -a $pooled_file -b $ind_file -wa -u > pooledRep/compInd/$ind_prefix.overlapped_pooledRep.bed
        done
done

pooled_files=($(ls pooledRep/*_pooledRep.bed))
for pooled_file in ${pooled_files[@]}; do
        prefix=$(basename $pooled_file | cut -d "_" -f 1-2)
        echo "$prefix"
        cat pooledRep/compInd/$prefix* | sort | uniq -c | awk 'BEGIN{OFS="\t"}{if($1 >=2) print $2,$3,$4,$5}' | sort-bed - > pooledRep/$prefix"_filtered.bed"
done
        
pooled_files=($(ls pooledRep/*_pooledRep.bed))
for pooled_file in ${pooled_files[@]}; do
        prefix=$(basename $pooled_file | cut -d "_" -f 1-2)
        before=$(wc -l pooledRep/$prefix"_pooledRep.bed" | awk '{print $1}') 
        after=$(wc -l pooledRep/$prefix"_filtered.bed" | awk '{print $1}')      
        echo -e "$prefix\t$before\t$after"
done > pooledRep/pooledRepPeakNum.txt

############### mergeCondition_pooledRep 
## MERGE ACROSS CONDITIONS
mkdir -p $atac_dir/peaks/mergeCondition_pooledRep
cd $atac_dir/peaks
bedtools merge -i <(cat pooledRep/*_filtered.bed | sort-bed -) -c 4 -o collapse -delim "|" > mergeCondition_pooledRep/all.filtered_withName.bed

## GENERAT SUMMARY FOR MERGE ALL
perl /mnt/isilon/sfgi/suc1/scripts/summarize_merged_bed.pl mergeCondition_pooledRep/all.filtered_withName.bed > mergeCondition_pooledRep/all.filtered_summary.txt

## RENAME MERGE ALL BED FILES
awk -v prefix=OCR -v i=1 'BEGIN{OFS="\t"}{name=prefix"_"i; i++; print $1,$2,$3,name}' mergeCondition_pooledRep/all.filtered_withName.bed > mergeCondition_pooledRep/all.filtered_newName.bed

## track
mkdir -p $atac_dir/peaks/track
files=($(ls pooledRep/*_filtered.bed))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "_" -f 1-2)
        if [[ $prefix == "CD4_unstim" ]]; then
                color="0,0,225"
        elif [[ $prefix == "CD4_8hr" ]]; then
                color="0,225,0"
        elif [[ $prefix == "CD4_24hr" ]]; then
                color="225,0,0"
        fi
        cp $file $atac_dir/peaks/track/$prefix"_filtered.track.bed"
        header=$(echo "track name=\"$prefix.atac\" description=\"$prefix filtered pooledPeak\" color=$color")
        echo $header
        sed -i "1i$header" $atac_dir/peaks/track/$prefix"_filtered.track.bed"
done

cp $atac_dir/peaks/mergeCondition_pooledRep/all.filtered_newName.bed $atac_dir/peaks/track/CD4_filtered.track.bed
header=$(echo "track name=\"CD4.atac\" description=\"all CD4 concensus ATAC\" color=0,0,0")
echo $header
sed -i "1i$header" $atac_dir/peaks/track/CD4_filtered.track.bed

cd  $atac_dir/peaks/track/
webserver="/var/www/html/ucsc/sfgi/suc1/wells/atacSeq/naiveT_PBMC"
scp -o PubkeyAuthentication=no -vvv *.bed suc1@reslncsfgweb01.research.chop.edu:$webserver



