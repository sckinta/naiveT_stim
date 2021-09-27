parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare" # change
experiment="naiveT_stim" # change
hicup_parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/hicup" # change
hicup_subdirs=("naiveT_unstimulated_ND517" "naiveT_unstimulated_TMP442" "naiveT_8hr_ND517" "naiveT_8hr_TMP442" "naiveT_24hr_ND517" "naiveT_24hr_TMP442") # change
juicer_subdirs=("naiveT_unstimulated_2reps" "naiveT_8hr_2reps" "naiveT_24hr_2reps") # change
juicer_parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer" # change
dir="$parent_dir/$experiment"


mkdir -p $dir

cd $dir/
mkdir -p $dir/data/
mkdir -p $dir/scripts/
mkdir -p $dir/plots/

###################### link files ####################
cd $dir/data/

resols=("40K" "10K" "4000" "2000" "1000")
for resol in ${resols[@]}; do
        mkdir -p $dir/data/$resol
        cd $dir/data/$resol
        for hicup_subdir in ${hicup_subdirs[@]}; do
                if [[ $resol == "1000" ]]; then
                        resol_old="1K"
                        file=$(ls $hicup_parent_dir/$hicup_subdir/matrix/*.$resol_old.cool)
                        newname=$(basename $file | sed "s/\.$resol_old\.cool/\.$resol\.cool/")
                        ln -s $file $newname
                else
                        file=$(ls $hicup_parent_dir/$hicup_subdir/matrix/*.$resol.cool)
                        ln -s $file
                fi
        done
        cd $dir/data
done

cd 

###################### HiCRep #######################
# 40K
files=($(ls $dir/data/40K/*.cool))
resol=40000
out_prefix="$dir/plots/scc_heatmap.40K"

echo "module load R/4.0.2
Rscript ~/hic/scripts/R/hicrep.R  $resol $out_prefix ${files[@]}
" > $dir/scripts/hicrep.40K.sh
cd $dir/scripts
qsub -cwd -l h_vmem=64G hicrep.40K.sh

# 10K
files=($(ls $dir/data/10K/*.cool))
resol=10000
out_prefix="$dir/plots/scc_heatmap.10K"

echo "module load R/4.0.2
Rscript ~/hic/scripts/R/hicrep.R  $resol $out_prefix ${files[@]}
" > $dir/scripts/hicrep.10K.sh
cd $dir/scripts
qsub -cwd -l h_vmem=128G hicrep.10K.sh

###################### multiHiCcompare loops #########################
# cooler dump by chr
resols=("4000" "2000" "1000")
for resol in ${resols[@]}; do
        # resol="10K"
        files=($(ls $dir/data/$resol/*.$resol.cool))
        chrs=(`seq 1 22` "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                mkdir -p $dir/data/$resol/$chr
                for file in ${files[@]}; do
                        prefix=$(basename $file | sed "s/\.cool//")
                        echo "source ~/.bashrc
cooler dump -r $chr -r2 $chr --join $file > $dir/data/$resol/$chr/$prefix.$chr.txt
" > $dir/scripts/cooler_dump.$prefix.$chr.sh
                        sed -i '1i#!/bin/bash' $dir/scripts/cooler_dump.$prefix.$chr.sh
                        cd $dir/scripts
                        sbatch --mem 8G -J $prefix.$chr cooler_dump.$prefix.$chr.sh
                        cd $dir
                done
        done
done

# check cooler dump successful
resols=("4000" "2000" "1000")
for resol in ${resols[@]}; do
        ls  $dir/data/$resol/*/*.txt | wc -l # 138
done

# run multiCompareHiC_byChr.R
resols=("1000" "2000" "4000")
resol="1000"
for resol in ${resols[@]}; do
        chrs=(`seq 1 22` "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                data_dir="$dir/data/$resol/$chr/"
                out_prefix="$dir/data/$resol/$chr/HiCcompare"
                numCores=1
                echo "module load R/4.0.2
        Rscript $dir/multiCompareHiC_byChr.R $data_dir $out_prefix $numCores
        " > $dir/data/$resol/$chr/multiCompareHiC.$resol.$chr.sh
                sed -i '1i#!/bin/bash' $dir/data/$resol/$chr/multiCompareHiC.$resol.$chr.sh
                cd $dir/data/$resol/$chr/
                # sbatch --mem 300G -t 24:00:00 multiCompareHiC.$resol.$chr.sh
                # cd $dir
        done
done

# check run multiCompareHiC_byChr.R successful
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        ls $dir/data/$resol/*/HiCcompare.data.Rdata | wc -l # 23
done

# resubmit failed run multiCompareHiC_byChr.R
resols=("1000" "2000" "4000")
for resol in ${resols[@]}; do
        chrs=(`seq 1 22` "X")
        for chr in ${chrs[@]}; do
                chr="chr"$chr
                if [[ ! -f $dir/data/$resol/$chr/HiCcompare.data.Rdata ]]; then
                        echo -e "$resol\t$chr"
                        # cd $dir/data/$resol/$chr/
                        # sbatch --mem 300G -t 24:00:00 multiCompareHiC.$resol.$chr.sh
                fi
        done
done

# recalculate fdr for multiCompareHiC
resols=("1000" "2000" "4000")
resols=("1000" "2000")
for resol in ${resols[@]}; do
        cd $dir/data/$resol
        echo "module load R/4.0.2
Rscript ~/hic/scripts/R/multiCompareHiC_recal_FDR.R $dir/data/$resol 0.05
" > multiCompareHiC_recal_FDR.sh
        sed -i '1i#!/bin/bash' multiCompareHiC_recal_FDR.sh
        sbatch --mem 300G -t 12:00:00 multiCompareHiC_recal_FDR.sh
done

# overlap with previous called consensus loop
multiCompareHiC_consensus_loop.R

######################## differential A/B compartments #################
mkdir -p $dir/data/40K/eigen
cd $dir/data/40K/eigen
conditions=("Acinar_2reps" "Alpha_2reps" "Beta_2reps")
for condition in ${conditions[@]}; do
        files=($(ls /mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/$condition/ABcompartments/*.eigen.40K.cis.vecs.tsv))
        for file in ${files[@]}; do
                ln -s $file
        done
done


        
