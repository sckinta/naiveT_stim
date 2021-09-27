atac_dir="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT"
cd $atac_dir
mkdir -p $atac_dir/annotation
cd $atac_dir/annotation
# promoter
mkdir -p $atac_dir/annotation/promoter
cd $atac_dir/annotation/promoter

ocr_bed="$atac_dir/peaks/mergeCondition_pooledRep/all.filtered_newName.bed"
gene_bed="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene.bed"
genetype_file="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt"
tss_bed="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.TSS.bed"
promoter_region="-1500,500"
out_dir="$atac_dir/annotation/promoter"

module load R/4.0.2
Rscript ~/scripts/R/annotate_prOCR.R $ocr_bed $gene_bed $genetype_file $tss_bed $promoter_region $out_dir