rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/Jurkat"
cd $rawdata_dir
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210122_A00303_0241_AHY5HLDRXX/Data/Intensities/BaseCalls/Jurkat*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        read=$(basename $file | cut -d "_" -f 1-3,5)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        ln -s $file $read.fastq.gz
        # echo $condition $read
done

rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT"
cd $rawdata_dir
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210122_A00303_0241_AHY5HLDRXX/Data/Intensities/BaseCalls/naive*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        read=$(basename $file | cut -d "_" -f 1-3,5)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        ln -s $file $read.fastq.gz
        # echo $condition $read
done


rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/Jurkat2"
cd $rawdata_dir
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210311_A00303_0259_AH3Y7FDRXY/Data/Intensities/BaseCalls/Jurkat*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        read=$(basename $file | cut -d "_" -f 1-3,5)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        ln -s $file $read.fastq.gz
        # echo $condition $read
done

rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT2"
cd $rawdata_dir
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210311_A00303_0259_AH3Y7FDRXY/Data/Intensities/BaseCalls/naiveT*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        read=$(basename $file | cut -d "_" -f 1-3,5)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        ln -s $file $read.fastq.gz
        # echo $condition $read
done

rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT3"
cd $rawdata_dir
files=($(ls /mnt/isilon/sfgi/novaseq_raw/210504_A00303_0277_BHC2M3DRXY/Data/Intensities/BaseCalls/CD4*))
for file in ${files[@]}; do
        condition=$(basename $file | cut -d "_" -f 1,2)
        read=$(basename $file | cut -d "_" -f 1-3,5)
        mkdir -p $rawdata_dir/$condition
        cd $rawdata_dir/$condition
        ln -s $file $read.fastq.gz
        # echo $condition $read
done

