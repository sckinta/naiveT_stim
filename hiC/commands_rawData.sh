rawdata_dir="/mnt/isilon/sfgi/rawData/wells/hiC/PMBC_naiveT_stim"

# link files
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210216_A00303_0248_BHTLKGDSXY/Data/Intensities/BaseCalls/naiveT*) $(ls /mnt/isilon/sfgi/novaseq_raw/210216_A00303_0249_AHTLKCDSXY/Data/Intensities/BaseCalls/naiveT*))

for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        sample=$(basename $file | cut -d "_" -f 1,2,3)
        read=$(basename $file | cut -d "_" -f 1,2,3,5)
        mkdir -p $rawdata_dir/$condition
        mkdir -p $rawdata_dir/$condition/$sample
        cd $rawdata_dir/$condition/$sample
        ln -s $file $read.fastq.gz
done


# count reads
cd $rawdata_dir
files=($(ls $rawdata_dir/*/*/*_R1.fastq.gz))
for file in ${files[@]}; do
        mydir=$(dirname $file)
        cd $mydir
        echo "zcat $file | wc -l | awk '{num=\$1/4; print num}' > $mydir/readCount.txt" > count_read.sh
        qsub -cwd count_read.sh
done


