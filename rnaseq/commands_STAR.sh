############################# Jurkat cells ###################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/Jurkat/STAR"
mkdir -p $dir
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/Jurkat"
condition_dirs=($(ls -d $rawdata_dir/*/))
for condition_dir in ${condition_dirs[@]}; do
        condition=$(basename $condition_dir | sed "s/\///")
        mkdir -p $dir/$condition
        cd $dir/$condition
        files=($(ls $condition_dir/*.fastq.gz))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "_" -f 1-3)
                mkdir -p $dir/$condition/$sample
                cd $dir/$condition/$sample
                ln -s $file
                # echo $condition $sample
        done
done

sample_dirs=($(ls -d $dir/*/*/))
for sample_dir in ${sample_dirs[@]}; do
        cd $sample_dir
        sample=$(basename $sample_dir | sed "s/\///")
        fastq1=$sample"_R1.fastq.gz"
        fastq2=$sample"_R2.fastq.gz"
        echo "source ~/.bashrc
cd $sample_dir
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesCommand zcat --readFilesIn $sample_dir/$fastq1 $sample_dir/$fastq2
samtools index Aligned.sortedByCoord.out.bam
sambamba flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
" > align_flagstat.sh
        qsub -cwd -l h_vmem=32G -o align.o -e align.e -N $sample align_flagstat.sh
done

cd $dir
files=($(ls */*/flagstat.txt))
for file in ${files[@]}; do
        condition=$(dirname $file | cut -d "/" -f 1)
        sample=$(dirname $file | cut -d "/" -f 2)
        alignCount=$(sed "1q;d" $file | awk '{print $1}')
        readPairCount=$(sed "7q;d" $file | awk '{print $1}')
        mappedCount=$(sed "5q;d" $file | awk '{print $1}')
        mappedPerc=$(sed "5q;d" $file | awk '{split($5,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        pairCount=$(sed "9q;d" $file | awk '{print $1}')
        pairPerc=$(sed "9q;d" $file | awk '{split($6,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        echo -e "$condition\t$sample\t$alignCount\t$readPairCount\t$mappedCount\t$mappedPerc\t$pairCount\t$pairPerc"
done > all.flagstat.txt
sed -i "1icondition\tsample\talignCount\treadPairCount\tmappedCount\tmappedPerc\tpairedCount\tpairedPerc" all.flagstat.txt


############################# PBMC_naiveT cells ###################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT/STAR"
mkdir -p $dir
cd $dir
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT"
condition_dirs=($(ls -d $rawdata_dir/*/))
for condition_dir in ${condition_dirs[@]}; do
        condition=$(basename $condition_dir | sed "s/\///")
        mkdir -p $dir/$condition
        cd $dir/$condition
        files=($(ls $condition_dir/*.fastq.gz))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "_" -f 1-3)
                mkdir -p $dir/$condition/$sample
                cd $dir/$condition/$sample
                ln -s $file
                # echo $condition $sample
        done
done

cd $dir
sample_dirs=($(ls -d $dir/*/*/))
for sample_dir in ${sample_dirs[@]}; do
        cd $sample_dir
        sample=$(basename $sample_dir | sed "s/\///")
        fastq1=$sample"_R1.fastq.gz"
        fastq2=$sample"_R2.fastq.gz"
        echo "source ~/.bashrc
cd $sample_dir
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesCommand zcat --readFilesIn $sample_dir/$fastq1 $sample_dir/$fastq2
samtools index Aligned.sortedByCoord.out.bam
sambamba flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
" > align_flagstat.sh
        qsub -cwd -l h_vmem=32G -o align.o -N $sample align_flagstat.sh
done

cd $dir
files=($(ls */*/flagstat.txt))
for file in ${files[@]}; do
        condition=$(dirname $file | cut -d "/" -f 1)
        sample=$(dirname $file | cut -d "/" -f 2)
        alignCount=$(sed "1q;d" $file | awk '{print $1}')
        readPairCount=$(sed "7q;d" $file | awk '{print $1}')
        mappedCount=$(sed "5q;d" $file | awk '{print $1}')
        mappedPerc=$(sed "5q;d" $file | awk '{split($5,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        pairCount=$(sed "9q;d" $file | awk '{print $1}')
        pairPerc=$(sed "9q;d" $file | awk '{split($6,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        echo -e "$condition\t$sample\t$alignCount\t$readPairCount\t$mappedCount\t$mappedPerc\t$pairCount\t$pairPerc"
done > all.flagstat.txt
sed -i "1icondition\tsample\talignCount\treadPairCount\tmappedCount\tmappedPerc\tpairedCount\tpairedPerc" all.flagstat.txt




############################# Jurkat cells2 ###################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/Jurkat2/STAR"
mkdir -p $dir
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/Jurkat2"
condition_dirs=($(ls -d $rawdata_dir/*/))
for condition_dir in ${condition_dirs[@]}; do
        condition=$(basename $condition_dir | sed "s/\///")
        mkdir -p $dir/$condition
        cd $dir/$condition
        files=($(ls $condition_dir/*.fastq.gz))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "_" -f 1-3)
                mkdir -p $dir/$condition/$sample
                cd $dir/$condition/$sample
                ln -s $file
                # echo $condition $sample
        done
done

sample_dirs=($(ls -d $dir/*/*/))
for sample_dir in ${sample_dirs[@]}; do
        cd $sample_dir
        sample=$(basename $sample_dir | sed "s/\///")
        fastq1=$sample"_R1.fastq.gz"
        fastq2=$sample"_R2.fastq.gz"
        echo "source ~/.bashrc
cd $sample_dir
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesCommand zcat --readFilesIn $sample_dir/$fastq1 $sample_dir/$fastq2
samtools index Aligned.sortedByCoord.out.bam
sambamba flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
" > align_flagstat.sh
        qsub -cwd -l h_vmem=32G -o align.o -e align.e -N $sample align_flagstat.sh
done

cd $dir
files=($(ls */*/flagstat.txt))
for file in ${files[@]}; do
        condition=$(dirname $file | cut -d "/" -f 1)
        sample=$(dirname $file | cut -d "/" -f 2)
        alignCount=$(sed "1q;d" $file | awk '{print $1}')
        readPairCount=$(sed "7q;d" $file | awk '{print $1}')
        mappedCount=$(sed "5q;d" $file | awk '{print $1}')
        mappedPerc=$(sed "5q;d" $file | awk '{split($5,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        pairCount=$(sed "9q;d" $file | awk '{print $1}')
        pairPerc=$(sed "9q;d" $file | awk '{split($6,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        echo -e "$condition\t$sample\t$alignCount\t$readPairCount\t$mappedCount\t$mappedPerc\t$pairCount\t$pairPerc"
done > all.flagstat.txt
sed -i "1icondition\tsample\talignCount\treadPairCount\tmappedCount\tmappedPerc\tpairedCount\tpairedPerc" all.flagstat.txt


############################# PBMC_naiveT cells ###################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT2/STAR"
mkdir -p $dir
cd $dir
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT2"
condition_dirs=($(ls -d $rawdata_dir/*/))
for condition_dir in ${condition_dirs[@]}; do
        condition=$(basename $condition_dir | sed "s/\///")
        mkdir -p $dir/$condition
        cd $dir/$condition
        files=($(ls $condition_dir/*.fastq.gz))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "_" -f 1-3)
                mkdir -p $dir/$condition/$sample
                cd $dir/$condition/$sample
                ln -s $file
                # echo $condition $sample
        done
done

cd $dir
sample_dirs=($(ls -d $dir/*/*/))
for sample_dir in ${sample_dirs[@]}; do
        cd $sample_dir
        sample=$(basename $sample_dir | sed "s/\///")
        fastq1=$sample"_R1.fastq.gz"
        fastq2=$sample"_R2.fastq.gz"
        echo "source ~/.bashrc
cd $sample_dir
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesCommand zcat --readFilesIn $sample_dir/$fastq1 $sample_dir/$fastq2
samtools index Aligned.sortedByCoord.out.bam
sambamba flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
" > align_flagstat.sh
        qsub -cwd -l h_vmem=32G -o align.o -N $sample align_flagstat.sh
done

cd $dir
files=($(ls */*/flagstat.txt))
for file in ${files[@]}; do
        condition=$(dirname $file | cut -d "/" -f 1)
        sample=$(dirname $file | cut -d "/" -f 2)
        alignCount=$(sed "1q;d" $file | awk '{print $1}')
        readPairCount=$(sed "7q;d" $file | awk '{print $1}')
        mappedCount=$(sed "5q;d" $file | awk '{print $1}')
        mappedPerc=$(sed "5q;d" $file | awk '{split($5,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        pairCount=$(sed "9q;d" $file | awk '{print $1}')
        pairPerc=$(sed "9q;d" $file | awk '{split($6,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        echo -e "$condition\t$sample\t$alignCount\t$readPairCount\t$mappedCount\t$mappedPerc\t$pairCount\t$pairPerc"
done > all.flagstat.txt
sed -i "1icondition\tsample\talignCount\treadPairCount\tmappedCount\tmappedPerc\tpairedCount\tpairedPerc" all.flagstat.txt


############################# PBMC_naiveT cells 3 ###################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3/STAR"
mkdir -p $dir
cd $dir
rawdata_dir="/mnt/isilon/sfgi/rawData/wells/rnaSeq/PBMC_naiveT3"
condition_dirs=($(ls -d $rawdata_dir/*/))
for condition_dir in ${condition_dirs[@]}; do
        condition=$(basename $condition_dir | sed "s/\///")
        mkdir -p $dir/$condition
        cd $dir/$condition
        files=($(ls $condition_dir/*.fastq.gz))
        for file in ${files[@]}; do
                sample=$(basename $file | cut -d "_" -f 1-3)
                mkdir -p $dir/$condition/$sample
                cd $dir/$condition/$sample
                ln -s $file
                # echo $condition $sample
        done
done

cd $dir
sample_dirs=($(ls -d $dir/*/*/))
for sample_dir in ${sample_dirs[@]}; do
        cd $sample_dir
        sample=$(basename $sample_dir | sed "s/\///")
        fastq1=$sample"_R1.fastq.gz"
        fastq2=$sample"_R2.fastq.gz"
        echo "source ~/.bashrc
cd $sample_dir
STAR --genomeDir /mnt/isilon/sfgi/indexes/star/hg19/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesCommand zcat --readFilesIn $sample_dir/$fastq1 $sample_dir/$fastq2
samtools index Aligned.sortedByCoord.out.bam
sambamba flagstat Aligned.sortedByCoord.out.bam > flagstat.txt
" > align_flagstat.sh
        sed -i '1i#!/bin/bash' align_flagstat.sh
        sbatch --mem 32G -t 24:00:00 -J $sample align_flagstat.sh
done

cd $dir
files=($(ls */*/flagstat.txt))
for file in ${files[@]}; do
        condition=$(dirname $file | cut -d "/" -f 1)
        sample=$(dirname $file | cut -d "/" -f 2)
        alignCount=$(sed "1q;d" $file | awk '{print $1}')
        readPairCount=$(sed "7q;d" $file | awk '{print $1}')
        mappedCount=$(sed "5q;d" $file | awk '{print $1}')
        mappedPerc=$(sed "5q;d" $file | awk '{split($5,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        pairCount=$(sed "9q;d" $file | awk '{print $1}')
        pairPerc=$(sed "9q;d" $file | awk '{split($6,a,":"); gsub(/\(/,"",a[1]); print a[1]}')
        echo -e "$condition\t$sample\t$alignCount\t$readPairCount\t$mappedCount\t$mappedPerc\t$pairCount\t$pairPerc"
done > all.flagstat.txt
sed -i "1icondition\tsample\talignCount\treadPairCount\tmappedCount\tmappedPerc\tpairedCount\tpairedPerc" all.flagstat.txt
