rawdata_dir="/mnt/isilon/sfgi/rawData/wells/captureC/PBMC_naiveT_unstim"
nova_dir="/mnt/isilon/sfgi/novaseq_raw/210318_A00303_0262_BHYCT3DMXX/Data/Intensities/BaseCalls"

files=($(ls $nova_dir/naiveT_unstimulated*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1-3)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        new_name_prefix=$(basename $file | cut -d "_" -f 1-3,5)
        ln -s $file $new_name_prefix.fastq.gz
done

# count read
cd $rawdata_dir
samples=($(ls -d */ | sed "s/\///"))
for sample in ${samples[@]}; do
        cd $rawdata_dir/$sample
        R1=$sample"_R1.fastq.gz"
        echo "zcat $R1 | wc -l  | awk '{num=\$1/4; print num}' > readCount.txt" > count_read.sh
        qsub -cwd count_read.sh
done


############### Jurkat ###################
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/captureC/Jurkat"
nova_dir="/mnt/isilon/sfgi/novaseq_raw/210401_A00303_0270_AH3KWWDSX2/Data/Intensities/BaseCalls"

files=($(ls $nova_dir/Jurkat*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1-2)
        sample=$(basename $file | cut -d "_" -f 1-3)
        mkdir -p $rawdata_dir/$condition
        mkdir -p $rawdata_dir/$condition/$sample
        cd $rawdata_dir/$condition/$sample
        new_name_prefix=$(basename $file | cut -d "_" -f 1-3,5)
        ln -s $file $new_name_prefix.fastq.gz
done

# count read
cd $rawdata_dir
mydirs=($(ls -d */* | sed "s/\/$//" | grep -v "Jurkat_unstim_rep3"))
for mydir in ${mydirs[@]}; do
        cd $rawdata_dir/$mydir
        sample=$(basename $(pwd))
        R1=$rawdata_dir/$mydir/$sample"_R1.fastq.gz"
        # echo -e '#!/bin/sh' > count_read.sh
        echo "zcat $R1 | wc -l  | awk '{num=\$1/4; print num}' > readCount.txt" > count_read.sh
        sed -i '1i#!/bin/bash' count_read.sh
        sbatch count_read.sh
done
