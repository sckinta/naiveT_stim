library(tidyverse)
if (!require(parseIbed)){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "parseIbed")
}
library(parseIbed)

if (!require(DFbedtools)){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "DFbedtools")
}
library(DFbedtools)

juicer_subdirs=("naiveT_unstimulated_2reps" "naiveT_8hr_2reps" "naiveT_24hr_2reps") # change
juicer_parent_dir="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer" # change

loop_lists="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.1500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.1kb.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.2000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.2500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.4000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps/mustache_ICE_loops/mustache_ICE_loops.5kb.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.1500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.1kb.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.2000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.2500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.4000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/mustache_ICE_loops/mustache_ICE_loops.5kb.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.1500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.1kb.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.2000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.2500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.4000.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.500.bedpe,/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps/mustache_ICE_loops/mustache_ICE_loops.5kb.bedpe"

tss_file="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.TSS.bed"
out_prefix="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/gene2OCR/Alpha_ArimaHIC_DpnII_Virtual.2reps"

loop_lists = strsplit(loop_lists,",")[[1]]

filter_by_distance=1
cis_only=1
prom_region=c(-1500,+500)

tss_bed=read_delim(
        tss_file, delim="\t", comment="#",
        col_names=c("tss_chr","tss_start","tss_end", "tss_anno","dot","strand")
) %>% 
mutate(tss_start=tss_start+1) %>% 
select(-dot)

## prom_bed
prom_bed = tss_bed %>% 
mutate(prom_start=ifelse(strand=="+", tss_start+prom_region[1], tss_start-prom_region[2]+1)) %>% 
mutate(prom_end=ifelse(strand=="+", tss_end+prom_region[2]-1, tss_end-prom_region[1])) %>% 
select(prom_chr=tss_chr, prom_start, prom_end, prom_anno=tss_anno) %>% 
separate_rows(prom_anno, sep="\\|")


read_bedpe <- function(bedpe){
        df = read_delim(
                bedpe, delim="\t", comment="#", 
                col_names=c("a_chr", "a_start", "a_end", "b_chr", "b_start", "b_end", "val1","val2")
        ) %>% 
        mutate(a_start=a_start+1) %>% 
        mutate(b_start=b_start+1) %>% 
        mutate(int_id=row_number())
        if (cis_only){
                df = df %>% filter(a_chr==b_chr)
        }
        
        baits = df %>% select(contains("a_")) %>% distinct()
        bait2prom = overlap_df(baits, prom_bed, minoverlap=1L)
        bait_anno = bind_cols(
                bait2prom$overlap_df1 %>% 
                select(bait_chr=a_chr, bait_start=a_start, bait_end=a_end),
                bait2prom$overlap_df2 %>% 
                select(prom_anno)
        )
        
        new_df = semi_join(
                df %>%
                dplyr::rename(bait_chr=a_chr, bait_start=a_start, bait_end=a_end, otherEnd_chr=b_chr, otherEnd_start=b_start, otherEnd_end=b_end), 
                bait_anno
        ) %>% 
        left_join(bait_anno)
        
        baits = df %>% select(contains("b_")) %>% distinct()
        bait2prom = overlap_df(baits, prom_bed, minoverlap=1L)
        bait_anno = bind_cols(
                bait2prom$overlap_df1 %>% 
                select(bait_chr=b_chr, bait_start=b_start, bait_end=b_end),
                bait2prom$overlap_df2 %>% 
                select(prom_anno)
        )
        
        new_df = bind_rows(
                new_df,
                semi_join(
                        df %>%
                        dplyr::rename(bait_chr=b_chr, bait_start=b_start, bait_end=b_end, otherEnd_chr=a_chr, otherEnd_start=a_start, otherEnd_end=a_end), 
                        bait_anno
                ) %>% 
                left_join(bait_anno)
        )
        new_df
}

anno_loops = lapply(
        loop_lists,
        read_bedpe
)