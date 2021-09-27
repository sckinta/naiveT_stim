parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer" # change
condition="naiveT_8hr_2reps" # change
dir="$parent_dir/$condition"
hicup_parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/hicup" # change
hicup_subdirs=("naiveT_8hr_ND517" "naiveT_8hr_TMP442") # change
cd $dir

mkdir -p $dir
cd $dir
mkdir -p $dir/scripts
# link cool files
mkdir -p $dir/cool

##### link cool files
cd $dir/cool
for hicup_subdir in ${hicup_subdirs[@]}; do
        files=($(ls $hicup_parent_dir/$hicup_subdir/matrix/*.cool))
        for file in ${files[@]};do
                ln -s $file
        done
done

###### merge cool
mkdir -p $dir/cool/merge

resolutions=("2.5M" "1M" "500K" "250K" "100K" "50K" "40K" "25K" "10K" "5K" "4000" "2500" "2000" "1500" "1K" "500")
resolutions=("1500" "500" "1K" "2000")
for resol in ${resolutions[@]}; do
        files=($(ls $dir/cool/*.$resol.cool))
        echo "conda activate py37
cooler merge $dir/cool/merge/$condition.$resol.cool ${files[@]}
cooler balance $dir/cool/merge/$condition.$resol.cool
" > $dir/scripts/cooler_merge.$resol.sh
        cd $dir/scripts
        qsub -cwd -l h_vmem=16G cooler_merge.$resol.sh
done

ln -s naiveT_8hr_2reps.1K.cool naiveT_8hr_2reps.1000.cool

####### A/B compartments - cooltools call-compartments
# https://github.com/open2c/cooltools/issues/116
# --reference-track /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.bed

# conda activate py37
# cooltools genome binnify /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes 40000 > /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins.40K.bed
# cooltools genome gc /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins.40K.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.fa > /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.bed
# cooltools genome genecov /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins.40K.bed hg19 > /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_genecov.40K.bed

mkdir -p $dir/ABcompartments
cd $dir/ABcompartments
files=($(ls $dir/cool/*.40K.cool) $(ls $dir/cool/merge/*.40K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
        echo "conda activate py37
cd $dir/ABcompartments
cooltools call-compartments --bigwig --reference-track /mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/hic_helper/hg19.bins_gc.40K.bed -o $dir/ABcompartments/$prefix.eigen.40K $file 
" > $dir/scripts/eigen.$prefix.sh
        cd $dir/scripts
        qsub -cwd eigen.$prefix.sh
        cd $dir/ABcompartments
done

# TAD boundaries - cooltools diamond-insulation
mkdir -p $dir/insulation_TADs
### The window size ({}) has to be a multiple of the bin size {}
files=($(ls $dir/cool/*.10K.cool) $(ls $dir/cool/merge/*.10K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
        echo "conda activate py37
cd $dir/insulation_TADs
cooltools diamond-insulation $file 500000 > $dir/insulation_TADs/$prefix.bin10K.window500K.tsv
" > $dir/scripts/insulation.$prefix.sh
        cd $dir/scripts
        qsub -cwd insulation.$prefix.sh
        cd $dir/insulation_TADs
done

# TADs - TADLib/hitad/output-DI (DI methods) - 10K
mkdir -p $dir/hitad_TADs
cd $dir/hitad_TADs

echo "res:10000" > $dir/hitad_TADs/meta_10k.txt
files=($(ls $dir/cool/*.10K.cool))
for i in `seq 1 ${#files[@]}`; do
        j=$((i-1))
        echo " rep$i:${files[$j]}"
done >>  $dir/hitad_TADs/meta_10k.txt

echo "conda activate py37
cd $dir/hitad_TADs
hitad -p 4 -O $dir/hitad_TADs/hiTADs_10K.bed -d $dir/hitad_TADs/meta_10k.txt --logFile $dir/hitad_TADs/hitad_10k.log
" > $dir/scripts/hitad_10k.sh
files=($(ls $dir/cool/*.10K.cool))
for file in ${files[@]}; do
        prefix=$(basename $file | cut -d "." -f 1)
echo "output-DI -O $dir/hitad_TADs/$prefix.10K.DI.bedGraph -p $file"
done >> $dir/scripts/hitad_10k.sh

cd $dir/scripts
qsub -cwd -pe smp 4 hitad_10k.sh

# # plot
# files=($(ls $dir/cool/*.40K.cool))
# for file in ${files[@]}; do
#         prefix=$(basename $file | cut -d "." -f 1)
# echo "conda activate py37
# cd $dir/hitad_TADs
# tad-plot -O $dir/hitad_TADs/$prefix.chr22.40K.png -p $file -T $dir/hitad_TADs/hiTADs_40K.bed -C chr22
# " > $dir/scripts/tad_plot.40k.$prefix.sh
# done


#################### interaction loops (mustache) ################################
mkdir -p $dir/mustache_ICE_loops
cd $dir/mustache_ICE_loops
resolutions=("500" "1K" "1500" "2000" "2500" "4000" "5K")
for resolution in ${resolutions[@]}; do
        res=$(echo $resolution | awk '{print tolower($0)}')
        if [[ $res =~ "k" || $res =~ "m" ]]; then
                res=$res"b"
        fi
        # echo "$res"
        mkdir -p $dir/mustache_ICE_loops/$res
        # cd $dir/mustache_ICE_loops/$res
        chrs=($(seq 1 22))
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                echo "conda activate mustache
mustache -f $dir/cool/merge/$condition.$resolution.cool -ch $chr -r $res -pt 0.1 -o $dir/mustache_ICE_loops/$res/loops.$chr.bedpe
" > $dir/scripts/mustache.$res.$chr.sh
                cd $dir/scripts
                qsub -cwd -l h_vmem=20G mustache.$res.$chr.sh
        done
done


cd $dir/mustache_ICE_loops
ls */*.bedpe | wc -l # 154
resolutions=($(ls -d */))
for resol in ${resolutions[@]}; do
        resol=$(echo $resol | sed "s/\///")
        files=($(ls $dir/mustache_ICE_loops/$resol/*.bedpe))
        for file in ${files[@]}; do
                sed 1d $file
        done > mustache_ICE_loops.$resol.bedpe
        # rm -rf $resol/
done

################################## fithic_loops #################################
mkdir -p $dir/fithic_loops
cd $dir/fithic_loops

# run fithic
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        mkdir -p $dir/fithic_loops/$resol
        cool_file=$dir/cool/merge/$condition.$resol.cool
        chrs=(`seq 1 22` "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                mkdir -p $dir/fithic_loops/$resol/$chr
                cd $dir/fithic_loops/$resol/$chr
                sbatch  -t 12:00:00 --mem 32G -J fithic ~/hic/scripts/bash/run_fithic.sh $chr $resol $cool_file $dir/fithic_loops/$resol/$chr
        done
        cd $dir/fithic_loops
done

ls */*/FitHiC.spline_pass2.res*.significances.txt.gz | wc -l # 69

# filter by fdr
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        # merge chr results
        echo "dir=$dir
resol=$resol
        " > $dir/fithic_loops/$resol/fithic_postprocess.sh
        echo 'chrs=($(seq 1 22) "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                if [ $chr = "chr1" ]; then
                        zcat $dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz
                else
                        zcat $dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz | sed 1d
                fi
        done | gzip > $dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.txt.gz
        ' >> $dir/fithic_loops/$resol/fithic_postprocess.sh

        # recal fdr
        infile="$dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.txt.gz"
        outdir="$dir/fithic_loops/$resol"
        echo "module load R/4.0.2
Rscript ~/hic/scripts/R/fithic_padjust.R $infile
        " >> $dir/fithic_loops/$resol/fithic_postprocess.sh
        
        # filter fdr < 1e-4
        filter_col="fdr"
        cutoff="1e-4"
        infile="$dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.significances.fdr.txt.gz"
        outfile="$dir/fithic_loops/$resol/FitHiC.spline_pass2.res$resol.sigInt.bedpe"
        echo "module load R/4.0.2
Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
        " >> $dir/fithic_loops/$resol/fithic_postprocess.sh

        sed -i '1i#!/bin/bash' $dir/fithic_loops/$resol/fithic_postprocess.sh
        cd $dir/fithic_loops/$resol
        sbatch -t 24:00:00 --mem 120G fithic_postprocess.sh
done

# hicACT pvalue adjust
resols=("1000" "2000" "4000")
chrs=($(seq 1 22) "X")
for resol in ${resols[@]}; do
        for chr in ${chrs[@]}; do
                # run hicACT
                chr="chr"$chr
                resolKB=$(echo $resol | sed "s/000//")
                infile="$dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz"
                outdir="$dir/fithic_loops/$resol/$chr"
                echo "module load R/4.0.2
Rscript ~/hic/scripts/R/hicACT.R $resolKB $infile $outdir
        " > $dir/fithic_loops/$resol/$chr/hicACT.sh
        
                # filter_fithic by ACT_adjusted < 1e-10
                filter_col="ACT_pvalue"
                cutoff="1e-10"
                infile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.txt.gz"
                outfile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe"
                echo "Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
" >> $dir/fithic_loops/$resol/$chr/hicACT.sh
        sed -i '1i#!/bin/bash' $dir/fithic_loops/$resol/$chr/hicACT.sh
        cd $dir/fithic_loops/$resol/$chr/
        sbatch --mem 64G -t 24:00:00 hicACT.sh
        cd $dir/fithic_loops
        done
done

# check failed HiACT
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        ls $dir/fithic_loops/$resol/*/ACT_adjusted.sigInt.bedpe| wc -l
done

resols=("1000" "2000" "4000")
chrs=($(seq 1 22) "X")
for resol in ${resols[@]}; do
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                if [ ! -f "$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe" ]; then
                        echo "$resol/$chr"
                fi
        done
done

# re-run failed HiACT
resols=("1000" "2000" "4000")
chrs=($(seq 1 22) "X")
for resol in ${resols[@]}; do
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                if [ ! -f "$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe" ]; then
                        # echo "$chr"
                        # run hicACT
                        resolKB=$(echo $resol | sed "s/000//")
                        infile="$dir/fithic_loops/$resol/$chr/FitHiC.spline_pass2.res$resol.significances.txt.gz"
                        outdir="$dir/fithic_loops/$resol/$chr"
                        echo "module load R/4.0.2
        Rscript ~/hic/scripts/R/hicACT.R $resolKB $infile $outdir
                " > $dir/fithic_loops/$resol/$chr/hicACT.sh

                        # filter_fithic by ACT_adjusted < 1e-10
                        filter_col="ACT_pvalue"
                        cutoff="1e-10"
                        infile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.txt.gz"
                        outfile="$dir/fithic_loops/$resol/$chr/ACT_adjusted.sigInt.bedpe"
                        echo "Rscript ~/hic/scripts/R/filter_fithic.R $filter_col $cutoff $resol $infile $outfile
        " >> $dir/fithic_loops/$resol/$chr/hicACT.sh
                sed -i '1i#!/bin/bash' $dir/fithic_loops/$resol/$chr/hicACT.sh
                cd $dir/fithic_loops/$resol/$chr/
                sbatch --mem 64G -t 24:00:00 hicACT.sh
                cd $dir/fithic_loops
                fi
        done
done



resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        cat $dir/fithic_loops/$resol/*/ACT_adjusted.sigInt.bedpe > $dir/fithic_loops/$resol/ACT_adjusted.res$resol.sigInt.bedpe
done

files=($(ls */*.bedpe))
for file in ${files[@]}; do
        wc -l $file
done

# 66558 1000/ACT_adjusted.res1000.sigInt.bedpe
# 78276 1000/FitHiC.spline_pass2.res1000.sigInt.bedpe
# 417286 2000/ACT_adjusted.res2000.sigInt.bedpe
# 548668 2000/FitHiC.spline_pass2.res2000.sigInt.bedpe
# 1285073 4000/ACT_adjusted.res4000.sigInt.bedpe
# 1756663 4000/FitHiC.spline_pass2.res4000.sigInt.bedpe

# further p-value cutoff
files=($(ls */ACT_adjusted.*.bedpe)) # ACT_pvalue < 1e-14
for file in ${files[@]}; do
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 < 1e-14) print $0}' $file | wc -l
done

files=($(ls */FitHiC.*.bedpe)) # FDR < 1e-6
for file in ${files[@]}; do
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 < 1e-6) print $0}' $file | wc -l
done
# 23702
# 239127
# 949103


############################# hic_loops ##############################
mkdir -p $dir/hic_loops
mkdir -p $dir/hic_loops/all
mkdir -p $dir/hic_loops/anno

### link loops
cd $dir/hic_loops/all
files=($(ls $dir/mustache_ICE_loops/*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 2)
        if [[ $resol == "1kb" ]]; then
                #echo "1000"
                ln -s $file $condition.mustache.pt01.res1000.bedpe
        elif [[ $resol == "2000" ]]; then
                #echo "2000"
                ln -s $file $condition.mustache.pt01.res2000.bedpe
        elif [[ $resol == "4000" ]]; then
                #echo "4000"
                ln -s $file $condition.mustache.pt01.res4000.bedpe
        fi
done

cd $dir/hic_loops/all
files=($(ls $dir/fithic_loops/*/ACT_adjusted.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 2)
        ln -s $file $condition.fithic_ACT.pv1e10.$resol.bedpe
done

cd $dir/hic_loops/all
files=($(ls $dir/fithic_loops/*/FitHiC.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 3)
        ln -s $file $condition.fithic_FDR.fdr1e4.$resol.bedpe 
done

### different cutoff for fithic_ACT and fithic_FDR
cd $dir/hic_loops/all
files=($(ls $dir/hic_loops/all/$condition.fithic_ACT.pv1e10.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 4)
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 <= 1e-14) print $0}' $file > $dir/hic_loops/all/$condition.fithic_ACT.pv1e14.$resol.bedpe
done

cd $dir/hic_loops/all
files=($(ls $dir/hic_loops/all/$condition.fithic_FDR.fdr1e4.*.bedpe))
for file in ${files[@]}; do
        resol=$(basename $file | cut -d "." -f 4)
        awk 'BEGIN{OFS="\t"; FS="\t"}{if($8 <= 1e-6) print $0}' $file > $dir/hic_loops/all/$condition.fithic_FDR.fdr1e6.$resol.bedpe
done

### generate annotated files - anno.bedpe, anno.gene2OCR.txt anno.summary.txt
mkdir -p $dir/hic_loops/scripts
# using concensus atac peak here
atac_peak_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/peaks/mergeCondition_pooledRep/all.filtered_newName.bed"

files=($(ls $dir/hic_loops/all/*.bedpe))
for file in ${files[@]}; do
        prefix=$(basename $file | sed "s/\.bedpe/\.anno/")
        outfile_prefix="$dir/hic_loops/anno/$prefix"
        echo "module load R/4.0.2
Rscript /mnt/isilon/sfgi/suc1/scripts/R/annotate_bedpe2geneOCR_hg19.R $file $atac_peak_file $outfile_prefix
" > $dir/hic_loops/scripts/$prefix.sh
        sed -i '1i#!/bin/bash' $dir/hic_loops/scripts/$prefix.sh
        cd $dir/hic_loops/scripts/
        sbatch --mem 8G $dir/hic_loops/scripts/$prefix.sh
done

### put summary together
cd $dir/hic_loops
files=($(ls $dir/hic_loops/anno/*.anno.summary.txt))
for i in `seq 1 ${#files[@]}`; do
        j=$((i-1))
        if [[ "$i" -eq "1" ]]; then
                cat ${files[$j]}
        else
                sed 1d ${files[$j]}
        fi
done > $dir/hic_loops/$condition.hic_loops.summary.txt

cd $dir/hic_loops
midfixes=($(ls all/*.bedpe | cut -d "/" -f 2 | cut -d "." -f 2-3 | sort -u))
for midfix in ${midfixes[@]}; do
        files=($(ls $dir/hic_loops/anno/*$midfix*.anno.gene2OCR.txt))
        geneOCR_pairN=$(cat ${files[@]} | sort -u | wc -l)
        gene_N=$(cat ${files[@]} | cut -f 3 | sort -u | wc -l)
        OCR_N=$(cat ${files[@]} | cut -f 1 | sort -u | wc -l)
        echo -e "$midfix\t$geneOCR_pairN\t$gene_N\t$OCR_N"
done > $dir/hic_loops/$condition.hic_loops_anno.noRes.summary.txt
sed -i "1i#prefix\tgeneOCR_pairN\tgene_N\tOCR_N" $dir/hic_loops/$condition.hic_loops_anno.noRes.summary.txt


        





