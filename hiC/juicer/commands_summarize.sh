juicer_subdirs=("naiveT_unstimulated_2reps" "naiveT_8hr_2reps" "naiveT_24hr_2reps") # change
juicer_parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer" # change
mkdir -p $juicer_parent_dir/naiveT_stim_summary
dir=$juicer_parent_dir/naiveT_stim_summary
cd $dir


## for anno
for condition in ${juicer_subdirs[@]}; do
        cat $juicer_parent_dir/$condition/hic_loops/anno/$condition.fithic_FDR.fdr1e6.*.anno.gene2OCR.txt $juicer_parent_dir/$condition/hic_loops/anno/$condition.mustache.pt01.*.anno.gene2OCR.txt
done | cut -f 1 | sort -u | wc -l
# 79930 total OCRs cross condition, resolution, caller

for condition in ${juicer_subdirs[@]}; do
        cat $juicer_parent_dir/$condition/hic_loops/anno/$condition.fithic_FDR.fdr1e6.*.anno.gene2OCR.txt $juicer_parent_dir/$condition/hic_loops/anno/$condition.mustache.pt01.*.anno.gene2OCR.txt
done | cut -f 3 | sort -u | wc -l
# 27544 total genes cross condition, resolution, caller

## for each resol
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        for condition in ${juicer_subdirs[@]}; do
                cat $juicer_parent_dir/$condition/hic_loops/all/$condition.fithic_FDR.fdr1e6.res$resol.bedpe $juicer_parent_dir/$condition/hic_loops/all/$condition.mustache.pt01.res$resol.bedpe | cut -f 1-6 
        done | sort -u | sort-bed - > $juicer_parent_dir/naiveT_stim_summary/naiveT_all.res$resol.bedpe
        wc -l $juicer_parent_dir/naiveT_stim_summary/naiveT_all.res$resol.bedpe
done

# 89779 /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_stim_summary/naiveT_all.res1000.bedpe
# 520227 /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_stim_summary/naiveT_all.res2000.bedpe
# 1612484 /mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_stim_summary/naiveT_all.res4000.bedpe

# resols=("1000" "2000" "4000")
# for resol in ${resols[@]}; do
#         for condition in ${juicer_subdirs[@]}; do
#                 cat $juicer_parent_dir/$condition/hic_loops/anno/$condition.fithic_FDR.fdr1e6.res$resol.anno.bedpe $juicer_parent_dir/$condition/hic_loops/anno/$condition.mustache.pt01.res$resol.anno.bedpe | cut -f 1-6,9-12 
#         done | sort -u > $juicer_parent_dir/naiveT_stim_summary/naiveT_all.res$resol.anno.bedpe
#         wc -l $juicer_parent_dir/naiveT_stim_summary/naiveT_all.res$resol.anno.bedpe
# done
