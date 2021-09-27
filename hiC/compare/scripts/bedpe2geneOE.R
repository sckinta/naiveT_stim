# R/3.6.0
library(readr)
library(dplyr)
library(tidyr)
library(parseIbed)
library(DFbedtools)

promoter_coords="-1500, 500"
search_promoter="both" # c("both","first","second")
tss_bed="/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.TSS.bed"
bedpe_file="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/juicer/Acinar_1/mustache_loops/ICE/mustache_loops_ICE.5kb.tsv"

bedpe2geneOE <- function(bedpe_file, search_promoter, tss_bed, promoter_coords){
        promoter_coords = as.integer(strsplit(promoter_coords,",")[[1]])
        if(grepl(".ibed", basename(bedpe_file))){
                bedpe = read_ibed_with_int_id(bedpe_file, parse_b2b="NULL", max_score=F) %>% 
                dplyr::rename(a_chr=bait_chr, a_start=bait_start, a_end=bait_end, b_chr=otherEnd_chr, b_start=otherEnd_start, b_end=otherEnd_end) %>% 
                mutate(bait_start=bait_start-1) %>% 
                mutate(otherEnd_start=otherEnd_start-1)
                
        }else if(grepl(".bedpe", basename(bedpe_file))){
                bedpe = vroom::vroom(bedpe_file, delim="\t", col_names = c("a_chr","a_start","a_end","b_chr","b_start","b_end"), comment="#")
        }else{
                bedpe = vroom::vroom(bedpe_file, delim="\t", col_names = c("a_chr","a_start","a_end","b_chr","b_start","b_end"), comment="#")
        }

        tss_df = vroom::vroom(tss_bed, delim="\t", col_names = c("tss_chr","tss_start","tss_pos","tss_anno","dot","strand"), comment="#") %>% 
        select(-tss_start, -dot)

        promoter_df = tss_df %>% 
        mutate(
                pro_start=ifelse(
                        strand=="+", tss_pos+promoter_coords[1], tss_pos-promoter_coords[2]
                )
        ) %>% 
        mutate(
                pro_end=ifelse(
                        strand=="+", tss_pos+promoter_coords[2], tss_pos-promoter_coords[1]
                )
        ) %>% dplyr::rename(pro_chr=tss_chr, pro_anno=tss_anno) %>% 
        select(pro_chr, pro_start, pro_end, pro_anno, strand)


        overlap_a = overlap_df(
                bedpe, promoter_df, 
                df1_chr_col="a_chr", df1_start_col="a_start", df1_end_col="a_end", df1_0base=F, 
                df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                minoverlap=1L
        )

        overlap_b = overlap_df(
                bedpe, promoter_df, 
                df1_chr_col="b_chr", df1_start_col="b_start", df1_end_col="b_end", df1_0base=F, 
                df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                minoverlap=1L
        )

        if (search_promoter=="both"){
                overlap_a = overlap_df(
                        bedpe, promoter_df, 
                        df1_chr_col="a_chr", df1_start_col="a_start", df1_end_col="a_end", df1_0base=F, 
                        df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                        minoverlap=1L
                )

                overlap_b = overlap_df(
                        bedpe, promoter_df, 
                        df1_chr_col="b_chr", df1_start_col="b_start", df1_end_col="b_end", df1_0base=F, 
                        df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                        minoverlap=1L
                )
        }else if(search_promoter=="first"){
                overlap_a = overlap_df(
                        bedpe, promoter_df, 
                        df1_chr_col="a_chr", df1_start_col="a_start", df1_end_col="a_end", df1_0base=F, 
                        df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                        minoverlap=1L
                )
        }else if(search_promoter=="second"){
                overlap_b = overlap_df(
                        bedpe, promoter_df, 
                        df1_chr_col="b_chr", df1_start_col="b_start", df1_end_col="b_end", df1_0base=F, 
                        df2_chr_col="pro_chr", df2_start_col="pro_start", df2_end_col="pro_end", df2_0base=F, 
                        minoverlap=1L
                )
        }


        if (search_promoter=="both"){
                pro2oe = bind_rows(
                        bind_cols(
                                overlap_a$overlap_df1,
                                overlap_a$overlap_df2,
                        ) %>% distinct(),
                        bind_cols(
                                overlap_b$overlap_df1 %>% 
                                dplyr::rename(b_chr=a_chr, b_start=a_start, b_end=a_end, a_chr=b_chr, a_start=b_start, a_end=b_end),
                                overlap_b$overlap_df2,
                        ) %>% distinct()
                ) %>% 
                dplyr::rename(oe_chr=b_chr, oe_start=b_start, oe_end=b_end) %>% 
                select(pro_anno, pro_chr, pro_start, pro_end, contains("oe_"), contains("a_")) %>% 
                distinct()
        }else if(search_promoter=="first"){
                pro2oe = bind_cols(
                        overlap_a$overlap_df1,
                        overlap_a$overlap_df2,
                ) %>% 
                dplyr::rename(oe_chr=b_chr, oe_start=b_start, oe_end=b_end) %>% 
                select(pro_anno, pro_chr, pro_start, pro_end, contains("oe_"), contains("a_")) %>% 
                distinct()
        }else if(search_promoter=="second"){
                pro2oe = bind_cols(
                        overlap_b$overlap_df1 %>% 
                        dplyr::rename(b_chr=a_chr, b_start=a_start, b_end=a_end, a_chr=b_chr, a_start=b_start, a_end=b_end),
                        overlap_b$overlap_df2,
                ) %>% dplyr::rename(oe_chr=b_chr, oe_start=b_start, oe_end=b_end) %>% 
                select(pro_anno, pro_chr, pro_start, pro_end, contains("oe_"), contains("a_")) %>% 
                distinct()
        }
        pro2oe
}


lapply(
        bedpe2geneOE,
        
)
