#R/4.0.2
library(tidyverse)
if (! require("parseIbed")){
        devtools::install_github("sckinta/myRpackages", ref = "master", subdir = "parseIbed")
}
if (! require("vroom")){
        install.packages("vroom")
}
library(parseIbed)
library(vroom)


args=commandArgs(trailingOnly=T)

bedpe_file=args[1]
atac_bed_file=args[2]
outfile_prefix=args[3]

# bedpe_file="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/hic_loops/all/naiveT_8hr_2reps.fithic_ACT.pv1e14.res1000.bedpe"
# atac_bed_file="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/peaks/mergeCondition_pooledRep/all.filtered_newName.bed"
# outfile_prefix="/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps/hic_loops/anno/naiveT_8hr_2reps.fithic_ACT.pv1e14.res1000.anno"

bedpe = vroom(
        bedpe_file, 
        col_names=c("chr_a","start_a","end_a","chr_b","start_b","end_b","val1","val2"),
        comment="#"
)

ocr_bed=vroom(
        atac_bed_file, 
        col_names=c("chr","start","end","id"),
        comment="#"
)

if (all(is.na(ocr_bed$id))) {
        ocr_bed = ocr_bed %>% 
        mutate(id=paste0("peak_",row_number()))
}


bedpe_geneOCR = annotate_bedpe2geneOCR(bedpe, ocr_bed, prom_bed)

bedpe_geneOCR = bedpe_geneOCR %>% 
mutate_at(vars(contains("start|end")), function(x){gsub(" ", "", format(x,scientific=F))})

write.table(
        bedpe_geneOCR,
        file=paste0(outfile_prefix,".bedpe"),
        sep="\t", row.names=F, col.names=F, quote=F
)

gene2ocr = bind_rows(
        bedpe_geneOCR %>% 
        select(ocr=ocr_a, anno=anno_b),
        bedpe_geneOCR %>% 
        select(ocr=ocr_b, anno=anno_a)
) %>% 
filter(!is.na(ocr), !is.na(anno)) %>% 
left_join(
        prom_bed %>% 
        select(anno=pro_anno, gene_id) %>% 
        distinct()
) %>% 
separate(anno,c("gene_name","transcript"), sep="\\+") %>% 
select(-transcript) %>% 
distinct()

write.table(
        gene2ocr,
        file=paste0(outfile_prefix,".gene2OCR.txt"),
        sep="\t", row.names=F, col.names=F, quote=F
)

### summarize
summary_out = tibble(
        file=bedpe_file,
        total_LoopN=bedpe %>% 
        select(contains("chr"),contains("start"), contains("end")) %>% 
        distinct() %>% 
        count() %>% pull(n),
        anno_LoopN = bedpe_geneOCR %>% 
        select(contains("chr"),contains("start"), contains("end")) %>% 
        distinct() %>% 
        count() %>% pull(n),
        geneOCR_pairN = nrow(gene2ocr),
        gene_N = gene2ocr %>% distinct(gene_id) %>% count() %>% pull(n),
        ocr_N = gene2ocr %>% distinct(ocr) %>% count() %>% pull(n)
)

write.table(
        summary_out,
        file=paste0(outfile_prefix,".summary.txt"),
        sep="\t", row.names=F, quote=F
)







