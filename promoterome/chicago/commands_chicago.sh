dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome/chicago" # change
chicago_run="naiveT_unstimulated_2reps" # change
atacSeq_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/Tonsil_Bcells/peaks/mergeRep/nB.filtered_newName.bed" # change (optional)
hicup_dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome" # change
hicup_names=("naiveT_unstimulated_ND365" "naiveT_unstimulated_ND523") # change

mkdir $dir/$chicago_run
mkdir $dir/$chicago_run/1frag
mkdir $dir/$chicago_run/4frag

######### FEATURES PREP
# MAKE SIMBOYLIC OF ATAC-SEQ CONSERVATIVE PEAK TO CURRENT DIR (optional)
cd $dir/$chicago_run
echo -e "DUMMY\t$atacSeq_file" > feature.list

#### 1FRAG CHICAGO
cd $dir/$chicago_run/1frag/
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_1frag/*_1frag.chinput)
        ln -s $file
done

# RUN 1FRAG CHICAGO
files=$(ls $dir/$chicago_run/1frag/*.chinput | tr "\n" "," | sed "s/,$//")

qsub -l m_mem_free=300G -l h_vmem=300G -cwd \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes

#### 4FRAG CHICAGO
cd $dir/$chicago_run/4frag
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_4frag/*_4frag.chinput)
        ln -s $file
done
# RUN 4FRAG CHICAGO
files=$(ls $dir/$chicago_run/4frag/*.chinput | tr "\n" "," | sed "s/,$//")

qsub -l m_mem_free=150G -l h_vmem=150G -cwd \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes




#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 4FRAG/
cd $dir/$chicago_run/4frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
qsub -cwd -l h_vmem=60G /mnt/isilon/sfgi/suc1/scripts/bash/exportIbed.sh

# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_4frag.txt

# # MAKE 4frag TRACKS
# perl ~/captureC/scripts/perl/makeWashUfromIbed.v3.pl chicagoRes/data/chicagoRes_gencodeV19.ibed > $chicago_run"_4frag"
# 
# # GENERATE GZIP FILE AND INDEX FOR WashU TRACK AND MOVE IT TO reslndbhitan03
# bgzip $chicago_run"_4frag"
# tabix -p bed $chicago_run"_4frag.gz"

# webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/EBV_Bcells/"  # change
# scp *.gz* suc1@reslncsfgweb01.research.chop.edu:$webserver_dir



#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 1FRAG/
cd $dir/$chicago_run/1frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
qsub -cwd -l h_vmem=32G /home/suc1/scripts/bash/exportIbed.sh
# /home/suc1/scripts/R/exportIbed.R


# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_1frag.txt

# # MAKE 1frag TRACKS
# perl ~/captureC/scripts/perl/makeWashUfromIbed.v3.pl chicagoRes/data/chicagoRes_gencodeV19.ibed > $chicago_run"_1frag"
# 
# # GENERATE GZIP FILE AND INDEX FOR WashU TRACK AND MOVE IT TO reslndbhitan03
# bgzip $chicago_run"_1frag"
# tabix -p bed $chicago_run"_1frag.gz"

# webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/EBV_Bcells/"  # change
# scp *.gz* suc1@reslncsfgweb01.research.chop.edu:$webserver_dir



############# merge
mkdir $dir/$chicago_run/merge
cd $dir/$chicago_run/merge
ln -s ../1frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_1frag.ibed
ln -s ../4frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_4frag.ibed

# # GENERATE 12Col ibed (takes 1-2 hr)
# bash ~/captureC/scripts/bash/merge2resIbed.sh chicagoRes_1frag.ibed  chicagoRes_4frag.ibed chicagoRes_merge.12col.ibed
# 
# 
# # FILTER AND CREATE 10 col regular ibed
# grep -v "discard" chicagoRes_merge.12col.ibed | cut -f 1-10  > chicagoRes_merge.ibed
module load R/4.0.2
Rscript ~/captureC/scripts/R/mergeIbedResForViz2.R -i chicagoRes_1frag.ibed,chicagoRes_4frag.ibed -o chicagoRes_merge.ibed

# CHICAGO SUMMARY ON merge ibed
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_merge_design3/chicago_both.baitmap chicagoRes_merge.ibed > summaryChicagoByInteraction_merge.txt

# # MAKE WashU TRACKs and upload to webserver
# perl ~/captureC/scripts/perl/makeWashUfromIbed.v3.pl chicagoRes_merge.ibed > $chicago_run"_merge"
# 
# bgzip $chicago_run"_merge"
# tabix -p bed $chicago_run"_merge.gz"
# 
# webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/EBV_Bcells/"  # change
# scp -o PubkeyAuthentication=no -vvv *.gz* suc1@reslncsfgweb01.research.chop.edu:$webserver_dir

# MAKE UCSC tracks
perl ~/captureC/scripts/perl/ibed2ucscBigInteract.pl chicagoRes_merge.ibed "no" "#7A67EE" | sort -k1,1 -k2,2n > chicagoRes_merge.bed

bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 chicagoRes_merge.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $chicago_run"_merge.bb"

rm chicagoRes_merge.bed

# webserver_dir="/var/www/html/ucsc/sfgi/suc1/grant/captureC/Promoterome/HepG2_new_3reps"  # change
# scp -o PubkeyAuthentication=no -vvv *.bb suc1@reslncsfgweb01.research.chop.edu:$webserver_dir


############ put summary files to one
cd $dir/$chicago_run
perl ~/captureC/scripts/perl/concatenateSummary.pl <(ls */*summaryChicagoByInteraction_*) > $chicago_run"_summaryChicagoByInteraction.txt"






################################## Jurkat 4hr ###################################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome/chicago" # change
chicago_run="Jurkat_4hr_3reps" # change
atacSeq_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/Tonsil_Bcells/peaks/mergeRep/nB.filtered_newName.bed" # change (optional)
hicup_dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome" # change
hicup_names=("Jurkat_4hr_rep1" "Jurkat_4hr_rep2" "Jurkat_4hr_rep3") # change

mkdir $dir/$chicago_run
mkdir $dir/$chicago_run/1frag
mkdir $dir/$chicago_run/4frag

######### FEATURES PREP
# MAKE SIMBOYLIC OF ATAC-SEQ CONSERVATIVE PEAK TO CURRENT DIR (optional)
cd $dir/$chicago_run
echo -e "DUMMY\t$atacSeq_file" > feature.list

#### 1FRAG CHICAGO
cd $dir/$chicago_run/1frag/
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_1frag/*_1frag.chinput)
        ln -s $file
done

# RUN 1FRAG CHICAGO
files=$(ls $dir/$chicago_run/1frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 80G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes

#### 4FRAG CHICAGO
cd $dir/$chicago_run/4frag
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_4frag/*_4frag.chinput)
        ln -s $file
done

# RUN 4FRAG CHICAGO
files=$(ls $dir/$chicago_run/4frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 60G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes




#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 4FRAG/
cd $dir/$chicago_run/4frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh

# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_4frag.txt


#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 1FRAG/
cd $dir/$chicago_run/1frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh
# /home/suc1/scripts/R/exportIbed.R


# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_1frag.txt


############# merge
mkdir $dir/$chicago_run/merge
cd $dir/$chicago_run/merge
ln -s ../1frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_1frag.ibed
ln -s ../4frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_4frag.ibed

# # GENERATE 12Col ibed (takes 1-2 hr)
# bash ~/captureC/scripts/bash/merge2resIbed.sh chicagoRes_1frag.ibed  chicagoRes_4frag.ibed chicagoRes_merge.12col.ibed
# 
# 
# # FILTER AND CREATE 10 col regular ibed
# grep -v "discard" chicagoRes_merge.12col.ibed | cut -f 1-10  > chicagoRes_merge.ibed
module load R/4.0.2
Rscript ~/captureC/scripts/R/mergeIbedResForViz2.R -i chicagoRes_1frag.ibed,chicagoRes_4frag.ibed -o chicagoRes_merge.ibed

# CHICAGO SUMMARY ON merge ibed
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_merge_design3/chicago_both.baitmap chicagoRes_merge.ibed > summaryChicagoByInteraction_merge.txt

# # MAKE WashU TRACKs and upload to webserver
# perl ~/captureC/scripts/perl/makeWashUfromIbed.v3.pl chicagoRes_merge.ibed > $chicago_run"_merge"
# 
# bgzip $chicago_run"_merge"
# tabix -p bed $chicago_run"_merge.gz"
# 
# webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/EBV_Bcells/"  # change
# scp -o PubkeyAuthentication=no -vvv *.gz* suc1@reslncsfgweb01.research.chop.edu:$webserver_dir

# MAKE UCSC tracks
perl ~/captureC/scripts/perl/ibed2ucscBigInteract.pl chicagoRes_merge.ibed "no" "#7A67EE" | sort -k1,1 -k2,2n > tmp.bed

bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 tmp.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $chicago_run"_merge.bb"

rm tmp.bed

webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/Jurkat"  # change
scp -o PubkeyAuthentication=no -vvv *.bb suc1@reslncsfgweb01.research.chop.edu:$webserver_dir


############ put summary files to one
cd $dir/$chicago_run
perl ~/captureC/scripts/perl/concatenateSummary.pl <(ls */*summaryChicagoByInteraction_*) > $chicago_run"_summaryChicagoByInteraction.txt"


################################## Jurkat 18hr ###################################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome/chicago" # change
chicago_run="Jurkat_18hr_3reps" # change
atacSeq_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/Tonsil_Bcells/peaks/mergeRep/nB.filtered_newName.bed" # change (optional)
hicup_dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome" # change
hicup_names=("Jurkat_18hr_rep1" "Jurkat_18hr_rep2" "Jurkat_18hr_rep3") # change

mkdir $dir/$chicago_run
mkdir $dir/$chicago_run/1frag
mkdir $dir/$chicago_run/4frag

######### FEATURES PREP
# MAKE SIMBOYLIC OF ATAC-SEQ CONSERVATIVE PEAK TO CURRENT DIR (optional)
cd $dir/$chicago_run
echo -e "DUMMY\t$atacSeq_file" > feature.list

#### 1FRAG CHICAGO
cd $dir/$chicago_run/1frag/
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_1frag/*_1frag.chinput)
        ln -s $file
done

# RUN 1FRAG CHICAGO
files=$(ls $dir/$chicago_run/1frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 80G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes

#### 4FRAG CHICAGO
cd $dir/$chicago_run/4frag
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_4frag/*_4frag.chinput)
        ln -s $file
done

# RUN 4FRAG CHICAGO
files=$(ls $dir/$chicago_run/4frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 60G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes




#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 4FRAG/
cd $dir/$chicago_run/4frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh

# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_4frag.txt


#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 1FRAG/
cd $dir/$chicago_run/1frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh
# /home/suc1/scripts/R/exportIbed.R


# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_1frag.txt


############# merge
mkdir $dir/$chicago_run/merge
cd $dir/$chicago_run/merge
ln -s ../1frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_1frag.ibed
ln -s ../4frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_4frag.ibed

# # GENERATE 12Col ibed (takes 1-2 hr)
# bash ~/captureC/scripts/bash/merge2resIbed.sh chicagoRes_1frag.ibed  chicagoRes_4frag.ibed chicagoRes_merge.12col.ibed
# 
# 
# # FILTER AND CREATE 10 col regular ibed
# grep -v "discard" chicagoRes_merge.12col.ibed | cut -f 1-10  > chicagoRes_merge.ibed
module load R/4.0.2
Rscript ~/captureC/scripts/R/mergeIbedResForViz2.R -i chicagoRes_1frag.ibed,chicagoRes_4frag.ibed -o chicagoRes_merge.ibed

# CHICAGO SUMMARY ON merge ibed
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_merge_design3/chicago_both.baitmap chicagoRes_merge.ibed > summaryChicagoByInteraction_merge.txt

# # MAKE WashU TRACKs and upload to webserver
# perl ~/captureC/scripts/perl/makeWashUfromIbed.v3.pl chicagoRes_merge.ibed > $chicago_run"_merge"
# 
# bgzip $chicago_run"_merge"
# tabix -p bed $chicago_run"_merge.gz"
# 
# webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/EBV_Bcells/"  # change
# scp -o PubkeyAuthentication=no -vvv *.gz* suc1@reslncsfgweb01.research.chop.edu:$webserver_dir

# MAKE UCSC tracks
perl ~/captureC/scripts/perl/ibed2ucscBigInteract.pl chicagoRes_merge.ibed "no" "#7A67EE" | sort -k1,1 -k2,2n > tmp.bed

bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 tmp.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $chicago_run"_merge.bb"

rm tmp.bed

webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/Jurkat"  # change
scp -o PubkeyAuthentication=no -vvv *.bb suc1@reslncsfgweb01.research.chop.edu:$webserver_dir


############ put summary files to one
cd $dir/$chicago_run
perl ~/captureC/scripts/perl/concatenateSummary.pl <(ls */*summaryChicagoByInteraction_*) > $chicago_run"_summaryChicagoByInteraction.txt"


################################## Jurkat unstim ###################################
dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome/chicago" # change
chicago_run="Jurkat_unstim_3reps" # change
atacSeq_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/Tonsil_Bcells/peaks/mergeRep/nB.filtered_newName.bed" # change (optional)
hicup_dir="/mnt/isilon/sfgi/suc1/analyses/wells/captureC/Promoterome" # change
hicup_names=("Jurkat_unstim_rep1" "Jurkat_unstim_rep2" "Jurkat_unstim_rep3") # change

mkdir $dir/$chicago_run
mkdir $dir/$chicago_run/1frag
mkdir $dir/$chicago_run/4frag

######### FEATURES PREP
# MAKE SIMBOYLIC OF ATAC-SEQ CONSERVATIVE PEAK TO CURRENT DIR (optional)
cd $dir/$chicago_run
echo -e "DUMMY\t$atacSeq_file" > feature.list

#### 1FRAG CHICAGO
cd $dir/$chicago_run/1frag/
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_1frag/*_1frag.chinput)
        ln -s $file
done

# RUN 1FRAG CHICAGO
files=$(ls $dir/$chicago_run/1frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 80G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_1frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes

#### 4FRAG CHICAGO
cd $dir/$chicago_run/4frag
# MAKE SIMBOYLIC LINK OF HICUP CHICAGO INPUT FILE TO CURRENT DIR
for hicup_name in ${hicup_names[@]}; do
        file=$(ls $hicup_dir/$hicup_name/hicup/*_4frag/*_4frag.chinput)
        ln -s $file
done

# RUN 4FRAG CHICAGO
files=$(ls $dir/$chicago_run/4frag/*.chinput | tr "\n" "," | sed "s/,$//")

sbatch --mem-per-cpu 60G -c 4 \
~/captureC/scripts/bash/runChicago.sh \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/settingsFile.txt \
/mnt/isilon/sfgi/manduchie/analyses/grant/captureC/PromoteromeDesign2/chicago_4frag_design2/ \
$dir/$chicago_run/feature.list \
$files chicagoRes


#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 4FRAG/
cd $dir/$chicago_run/4frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh

# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_4frag_design3/chicago_4frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_4frag.txt


#### GENERATE CHICAGO SUMMARY AND WASHU TRACKS IN 1FRAG/
cd $dir/$chicago_run/1frag/

# GENERATE IBED FILE IN R (exportIbed.sh used exportIbed.R which both saved in my home scripts)
sbatch --mem 32G ~/captureC/scripts/bash/exportIbed.sh
# /home/suc1/scripts/R/exportIbed.R


# CHANGE IBED ANNOTATION TO GENCODE_V19
perl ~/captureC/scripts/perl/reannotateIbed.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes.ibed > chicagoRes/data/chicagoRes_gencodeV19.ibed

# GENERATE CHICAGO SUMMARY USING NEW scripts
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_1frag_design3/chicago_1frag.baitmap \
chicagoRes/data/chicagoRes_gencodeV19.ibed > summaryChicagoByInteraction_1frag.txt


############# merge
mkdir $dir/$chicago_run/merge
cd $dir/$chicago_run/merge
ln -s ../1frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_1frag.ibed
ln -s ../4frag/chicagoRes/data/chicagoRes_gencodeV19.ibed chicagoRes_4frag.ibed

# # GENERATE merge ibed
module load R/4.0.2
Rscript ~/captureC/scripts/R/mergeIbedResForViz2.R -i chicagoRes_1frag.ibed,chicagoRes_4frag.ibed -o chicagoRes_merge.ibed

# CHICAGO SUMMARY ON merge ibed
perl ~/captureC/scripts/perl/summaryChicagoByInteraction.pl \
/mnt/isilon/sfgi/suc1/analyses/grant/captureC/PromoteromeDesign_gencodeV19/chicago_merge_design3/chicago_both.baitmap chicagoRes_merge.ibed > summaryChicagoByInteraction_merge.txt

# MAKE UCSC tracks
perl ~/captureC/scripts/perl/ibed2ucscBigInteract.pl chicagoRes_merge.ibed "no" "#7A67EE" | sort -k1,1 -k2,2n > tmp.bed

bedToBigBed -as=/mnt/isilon/sfgi/referenceSequences/interact.as -type=bed5+13 tmp.bed /mnt/isilon/sfgi/referenceSequences/hg19/hg19.chrom.sizes $chicago_run"_merge.bb"

rm tmp.bed

webserver_dir="/var/www/html/ucsc/sfgi/suc1/wells/captureC/Promoterome/Jurkat"  # change
scp -o PubkeyAuthentication=no -vvv *.bb suc1@reslncsfgweb01.research.chop.edu:$webserver_dir


############ put summary files to one
cd $dir/$chicago_run
perl ~/captureC/scripts/perl/concatenateSummary.pl <(ls */*summaryChicagoByInteraction_*) > $chicago_run"_summaryChicagoByInteraction.txt"

