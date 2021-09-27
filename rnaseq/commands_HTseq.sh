###################### Jurkat cells ##########################
# NEBNext kit non-stranded -- no
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/Jurkat"
mkdir -p $rnaseq_dir/HTseq
cd $rnaseq_dir/HTseq
mkdir -p $rnaseq_dir/HTseq/scripts
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.no.txt"
        job_file=$sample"_htseq_no.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s no -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        qsub -cwd -l h_vmem=32G $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done

perl ~/scripts/mergeReadCountFile.pl <(ls *_HTSeq.no.txt) > htseq_counts.txt

###################### PBMC_naiveT cells ##########################
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.no.txt"
        job_file=$sample"_htseq_no.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s no -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        qsub -cwd -l h_vmem=32G $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done


perl ~/scripts/mergeReadCountFile.pl <(ls *_HTSeq.no.txt) > htseq_counts.txt


# new batch kit: NEBNext Ultra II Directional RNA Library Prep Kit for Illumina (Cat# E7760S)


###################### Jurkat cells2 ##########################
# NEBNext kit stranded -- reverse
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/Jurkat2"
mkdir -p $rnaseq_dir/HTseq
cd $rnaseq_dir/HTseq
mkdir -p $rnaseq_dir/HTseq/scripts
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.txt"
        job_file=$sample"_htseq.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s reverse -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        qsub -cwd -l h_vmem=32G $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done

perl ~/scripts/mergeReadCountFile.pl <(ls *_HTSeq.txt) > htseq_counts.txt

###################### PBMC_naiveT cells 2##########################
# NEBNext kit stranded -- reverse
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT2"
mkdir -p $rnaseq_dir/HTseq
cd $rnaseq_dir/HTseq
mkdir -p $rnaseq_dir/HTseq/scripts
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.txt"
        job_file=$sample"_htseq.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s reverse -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        qsub -cwd -l h_vmem=32G $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done

perl ~/scripts/mergeReadCountFile.pl <(ls *_HTSeq.txt) > htseq_counts.txt



rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT2"
mkdir -p $rnaseq_dir/HTseq
cd $rnaseq_dir/HTseq
mkdir -p $rnaseq_dir/HTseq/scripts
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.no.txt"
        job_file=$sample"_htseq.no.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s no -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.srt.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        qsub -cwd -l h_vmem=32G $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done


# NEBNext Ultra II Directional RNA Library Prep Kit for Illumina (Cat# E7760S). - stranded -- reverse
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3"
mkdir -p $rnaseq_dir/HTseq
cd $rnaseq_dir/HTseq
mkdir -p $rnaseq_dir/HTseq/scripts
bam_files=($(ls $rnaseq_dir/STAR/*/*/Aligned.sortedByCoord.out.bam))
for bam_file in ${bam_files[@]}; do
        sample=$(dirname $bam_file | rev | cut -d '/' -f 1 | rev)
        # echo $sample
        output=$sample"_HTSeq.txt"
        job_file=$sample"_htseq.sh"
        echo "module load python/2.7
htseq-count -f bam -r pos -s reverse -t exon -m union $bam_file  /mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencode.v19.annotation.gtf > $rnaseq_dir/HTseq/$output
" > $rnaseq_dir/HTseq/scripts/$job_file
        sed -i '1i#!/bin/bash' $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq/scripts/
        sbatch --mem 32G -t 24:00:00 $rnaseq_dir/HTseq/scripts/$job_file
        cd $rnaseq_dir/HTseq
done

perl ~/scripts/mergeReadCountFile.pl <(ls *_HTSeq.txt) > htseq_counts.txt