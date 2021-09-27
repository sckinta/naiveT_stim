atac_dir="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT"
cd $atac_dir
## ATAC-seq encode pipeline was performed by Matt P. for each replicate
## consistent individuals are ND517, ND578, TMP442

## link bam_filt
conditions=("CD4_unstim" "CD4_8hr" "CD4_24hr")
for condition in ${conditions[@]}; do
        mkdir -p $atac_dir/$condition
        cd $atac_dir/$condition
        condition_name=$(echo $condition | sed "s/CD4_//")
        samples=("ND517" "ND578" "TMP442")
        for sample in ${samples[@]}; do
                sample_name=$sample"_"$condition_name
                # ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample_name/out/
                bam_filt=$(ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample_name/out/align/rep1/*.PE2SE.nodup.bam)
                ln -s $bam_filt $condition"_"$sample.nodup.bam
                ln -s $bam_filt.bai $condition"_"$sample.nodup.bam.bai
        done
done

## prepare encode
conditions=("CD4_unstim" "CD4_8hr" "CD4_24hr")
for condition in ${conditions[@]}; do
        cd $atac_dir/$condition
        echo "use_system = slurm
q = all.q
mem_bwt2 = 60G
mem_dedup = 16G
mem_macs2 = 20G
mem_ataqc = 200G
trimmed_fastq = true
pe = true
nth = 8" > atac.config
        echo "title = $condition" >> atac.config
        echo "enable_idr = true
species = hg19
species_file = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/bds_atac_species.conf" >> atac.config
        samples=($(ls *.bam | cut -d "." -f 1))
        for i in `seq 1 ${#samples[@]}`; do
                j=$((i-1))
                echo "filt_bam$i"" = ${samples[$j]}.nodup.bam"
        done >> atac.config
        echo "tss_enrich = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz
dnase = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz
prom = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz
enh = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz
reg2map = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/dnase_avgs_reg2map_p10_merged_named.pvals.gz
roadmap_meta = /mnt/isilon/sfgi/atac_dnase_pipelines_genome_data/hg19/ataqc/eid_to_mnemonic.txt
"  >> atac.config
        # cd $atac_dir/
done

# run encode
conditions=("CD4_unstim" "CD4_8hr" "CD4_24hr")
for condition in ${conditions[@]}; do
        cd $atac_dir/$condition
        echo "bds -c /mnt/isilon/sfgi/programs/.bds/bds.config /mnt/isilon/sfgi/programs/atac_dnase_pipelines/atac.bds $atac_dir/$condition/atac.config >> bds.log 2>&1" > atac_slurm.sh
        sed -i '1i#!/bin/bash' atac_slurm.sh
        sbatch --mem 16G -t 72:00:00 atac_slurm.sh
done


# summarize encode qc
mkdir -p $atac_dir/qc
cd $atac_dir/qc
samples=("ND517" "ND578" "TMP442")
for sample in ${samples[@]}; do
        ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample*/out/qc/rep*/*.PE2SE_qc.txt
done > file_list.txt
perl ~/functionalGenomics/scripts/perl/extract_atac_qc.pl file_list.txt > $atac_dir/qc/QC_summary.txt
rm file_list.txt

mkdir -p $atac_dir/qc/insert_hist
cd $atac_dir/qc/insert_hist
for sample in ${samples[@]}; do
        files=($(ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample*/out/qc/rep*/*.inserts.hist_graph.pdf))
        for file in ${files[@]}; do
                ln -s $file
        done
done

mkdir -p $atac_dir/qc/xcor_plots
cd $atac_dir/qc/xcor_plots
for sample in ${samples[@]}; do
        files=($(ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample*/out/qc/rep*/*.PE2SE.nodup.tn5.no_chrM.25M.R1.cc.plot.png))
        for file in ${files[@]}; do
                ln -s $file
        done
done

samples=("ND517" "ND578" "TMP442")
for sample in ${samples[@]}; do
        files=($(ls /mnt/isilon/sfgi/pahlm/analyses/wells/atacSeq/Zack_Jurkat_Nov2020/$sample*/out/qc/rep*/*.PE2SE_qc.txt))
        for file in ${files[@]}; do
                FRinNFR=$(grep "Fraction of reads in NFR" $file | cut -f 2)
                FRinPeak=$(grep "Fraction of reads in called peak regions" $file | cut -f 3)
                prefix=$(basename $file | cut -d "." -f 1)
                echo -e "$prefix\t$FRinNFR\t$FRinPeak"
        done
done > $atac_dir/qc/QC_summary_additional.txt
sed -i '1iprefix\tFRinNFR\tFRinPeak' $atac_dir/qc/QC_summary_additional.txt


################################ bigwig ##################################
mkdir -p $atac_dir/bigwig
mkdir -p $atac_dir/bigwig/scripts

#### individuals bams
bams=($(ls $atac_dir/*/*.bam))
for bam in ${bams[@]}; do
        prefix=$(basename $bam | cut -d "." -f 1)
        echo "source ~/.bashrc
conda activate deeptools
bamCoverage --normalizeUsing RPGC --effectiveGenomeSize 2864785220 --binSize 10 --extendReads -b $bam -o $atac_dir/bigwig/$prefix.bw
" > $atac_dir/bigwig/scripts/bam2bw.$prefix.sh
        sed -i '1i#!/bin/bash' $atac_dir/bigwig/scripts/bam2bw.$prefix.sh
        cd $atac_dir/bigwig/scripts
        sbatch --mem 16G -t 12:00:00 bam2bw.$prefix.sh
done

cd $atac_dir/bigwig
webserver="/var/www/html/ucsc/sfgi/suc1/wells/atacSeq/naiveT_PBMC"
scp -o PubkeyAuthentication=no -vvv *.bw suc1@reslncsfgweb01.research.chop.edu:$webserver

files=($(ls *.bw))
for file in ${files[@]}; do
        prefix=$(echo $file | sed "s/\.bw//")
        echo -e "ATAC-seq bw\t$prefix\t/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/bigwig/$prefix.bw\t/var/www/html/ucsc/sfgi/suc1/wells/atacSeq/naiveT_PBMC/$prefix.bw\thttp://159.14.198.88/ucsc/sfgi/suc1/wells/atacSeq/naiveT_PBMC/$prefix.bw\ttrack type=bigWig name=\"$prefix\" description=\"$prefix\" color=165,183,49 visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off bigDataUrl=http://159.14.198.88/ucsc/sfgi/suc1/wells/atacSeq/naiveT_PBMC/$prefix.bw"
done > copy2store.txt

#### merge bams
conditions=("CD4_unstim" "CD4_8hr" "CD4_24hr")
for condition in ${conditions[@]}; do
        bams=($(ls $atac_dir/$condition/*.bam))
        echo "source ~/.bashrc
samtools merge $atac_dir/bigwig/$condition.merge.bam ${bams[@]}
samtools sort -o $atac_dir/bigwig/$condition.merge.srt.bam $atac_dir/bigwig/$condition.merge.bam
samtools index $atac_dir/bigwig/$condition.merge.srt.bam
conda activate deeptools
bamCoverage --normalizeUsing RPGC --effectiveGenomeSize 2864785220 --binSize 10 --extendReads -b $atac_dir/bigwig/$condition.merge.srt.bam -o $atac_dir/bigwig/$condition.merge.bw
" > $atac_dir/bigwig/scripts/merge_bam2bw.$condition.sh
        sed -i '1i#!/bin/bash' $atac_dir/bigwig/scripts/merge_bam2bw.$condition.sh
        cd $atac_dir/bigwig/scripts
        sbatch --mem 32G -t 12:00:00 $atac_dir/bigwig/scripts/merge_bam2bw.$condition.sh
done
        



