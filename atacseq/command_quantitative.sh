atac_dir="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT"
cd $atac_dir

mkdir -p $atac_dir/quantitative
mkdir -p $atac_dir/quantitative/scripts/
cd $atac_dir/quantitative

## designFile
conditions=("CD4_unstim" "CD4_8hr" "CD4_24hr")
for condition in ${conditions[@]}; do
        files=($(ls $atac_dir/$condition/*.bam))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "." -f 1)
                echo -e "$sample\t$condition\t$file"
        done
done > designFile.txt
sed -i "1isample\tcondition\tfilePath" designFile.txt


## normalize
species="hg19"
design_file="$atac_dir/quantitative/designFile.txt"
peak_bed="$atac_dir/peaks/mergeCondition_pooledRep/all.filtered_newName.bed"
out_prefix="$atac_dir/quantitative/all.filtered_newName"
echo "module load R/4.0.2
Rscript ~/functionalGenomics/scripts/R/csaw_normalize.R $species $design_file $peak_bed $out_prefix
" > $atac_dir/quantitative/scripts/csaw_normalize.sh
sed -i '1i#!/bin/bash' $atac_dir/quantitative/scripts/csaw_normalize.sh
cd $atac_dir/quantitative/scripts
sbatch --mem 32G -t 12:00:00 csaw_normalize.sh




