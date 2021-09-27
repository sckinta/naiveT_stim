# fragInt2geneOCR.R
# R/3.3.2

library(tidyverse)
if (!require(parseIbed)){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "parseIbed")
}
library(parseIbed)

if (!require(DFbedtools)){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "DFbedtools")
}
library(DFbedtools)

args=commandArgs(trailingOnly=T)

ocr_file=args[1]
loop_lists=args[2]
tss_file=args[3]
out_prefix=args[4]

# ocr_file="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/peaks/Alpha.repro.bed"
# loop_lists="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/loops/Alpha_ArimaHIC_DpnII_Virtual.2reps.1frag.ibed,/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/loops/Alpha_ArimaHIC_DpnII_Virtual.2reps.4frag.ibed"
# tss_file="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.TSS.bed"
# out_prefix="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/gene2OCR/Alpha_ArimaHIC_DpnII_Virtual.2reps"

loop_lists = strsplit(loop_lists,",")[[1]]

### parameters
filter_by_distance=1
cis_only=1
prom_region=c(-1500,+500)

## read ocr_file
ocr_bed=read_delim(
        ocr_file, delim="\t", 
        col_names=c("ocr_chr","ocr_start","ocr_end", "ocr_id")
) %>% 
mutate(ocr_start=ocr_start+1)


## read tss_file
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

## function read_ibed
# ibed=loop_lists[1]
read_ibed <- function(ibed){
        df = read_ibed_with_int_id(
                ibed, parse_b2b="only bi-direction", max_score=F
        )
        if (cis_only){
                df = df %>% filter(bait_chr==otherEnd_chr)
        }
        if (filter_by_distance){
                baits = df %>% select(contains("bait_")) %>% distinct()
                bait2prom = overlap_df(baits, prom_bed, minoverlap=1L)
                bait_anno = bind_cols(
                        bait2prom$overlap_df1 %>% 
                        select(bait_chr, bait_start, bait_end),
                        bait2prom$overlap_df2 %>% 
                        select(prom_anno)
                )
                
                new_df = semi_join(df, bait_anno) %>% 
                select(-bait_name) %>% 
                left_join(bait_anno)
        }else{
                new_df = df %>% separate_rows(bait_name, sep="\\|") %>% 
                dplyr::rename(prom_anno=bait_name)
        }
        new_df %>% select(-otherEnd_name)
}


## function read_bedpe
# bedpe="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/virtualComp/arima_benchmark/loops/HepG2_ArimaHIC_loops.2reps.1000bp.bedpe"
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


## bait2gene
bait2gene = lapply(
        loop_lists,
        function(file){
                if(grepl("\\.ibed$",basename(file))){
                        read_ibed(file) %>% 
                        mutate(file=file)
                }else if(grepl("\\.bedpe$",basename(file))){
                        read_bedpe(file) %>% 
                        mutate(file)
                }
        }
)

bait2gene = do.call("bind_rows",bait2gene)

## otherEnd2OCR

otherEnds = bait2gene %>% 
select(contains("otherEnd_")) %>% 
distinct()

otherEnds2OCR=overlap_df(otherEnds, ocr_bed, minoverlap=1L)

otherEnds2OCR = bind_cols(
        otherEnds2OCR$overlap_df1,
        otherEnds2OCR$overlap_df2
)

## fragInt2geneOCR

fragInt2geneOCR = semi_join(bait2gene, otherEnds2OCR) %>% 
left_join(otherEnds2OCR)

save(fragInt2geneOCR, file=paste0(out_prefix, ".Rdata"))
