# R/3.3.2
library(readr)
library(dplyr)
library(tidyr)

########################### PBMC_naiveT ##############################
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3"

########### read count column and row filter ##################
geneCount=read_delim(file.path(rnaseq_dir,"HTseq","htseq_counts.txt"), delim="\t")
geneCount=geneCount %>% filter(!grepl("__",gene))
colnames(geneCount)=gsub("_HTSeq","",colnames(geneCount))

readCount_summary = tibble(
        sample=colnames(geneCount %>% select(-gene)),
        gene_align_count=apply(geneCount %>% select(-gene), 2, sum)
)

# read chrM and rRNA 
gene_id_rm = bind_rows(
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencode.v19.chrM_gene_id.txt", 
                delim="\t", col_names=c("gene_id","type")
        ),
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencode.v19.rRNA_gene_id.txt", 
                delim="\t", col_names=c("gene_id","type")
        )
) %>% distinct(gene_id)

# remove chrM and rRNA
geneCount = geneCount %>% filter(!gene %in% gene_id_rm$gene_id)
readCount_summary = left_join(
        readCount_summary,
        tibble(
                sample=colnames(geneCount %>% select(-gene)),
                gene_rm_chrM_rRNA_count=apply(geneCount %>% select(-gene), 2, sum)
        )
)

# just keep non-pseudogene and add gene_name
gene_type = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#")

geneCount = geneCount %>% 
dplyr::rename(gene_id=gene) %>% 
left_join(
        gene_type %>% 
        select(gene_id, gene_name, gene_type, biotype) %>% 
        distinct()
) %>% filter(!gene_type %in% c("pseudogene"))

readCount_summary = left_join(
        readCount_summary,
        tibble(
                sample=colnames(geneCount %>% select(-gene_name, -gene_id,-gene_type, -biotype)),
                gene_nonpseudo_count=apply(geneCount %>% select(-gene_name, -gene_id,-gene_type, -biotype), 2, sum)
        )
)

# remove genes using TPM < 1
gene_len=read_delim("/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencode.v19.gene_length.txt",delim="\t", 
col_names=c("gene_id","gene_name","gene_len"))

geneCount_melt = left_join(
        geneCount %>% 
        gather(key="sample",value="rawCount", -gene_id, -gene_name, -gene_type, -biotype),
        gene_len
)
geneCount_melt = geneCount_melt %>% mutate(RPK=rawCount/gene_len*1000)

scale_factor=geneCount_melt %>% 
group_by(sample) %>% 
summarise(total_RPK=sum(RPK)) %>% 
ungroup() %>% 
mutate(scale_factor=total_RPK/1000000)

tpm=left_join(geneCount_melt, scale_factor) %>% 
mutate(tpm=RPK/scale_factor) %>% 
select(gene_id,gene_name,gene_type, sample,tpm)

tpm = tpm %>% 
separate(sample, c("cell","time","ind"), sep="_", remove=F) %>% 
unite("condition",cell, time, sep="_")

meanTPM=tpm %>% 
group_by(gene_id, gene_name, condition) %>% 
summarise(tpm=mean(tpm)) %>% 
ungroup()

exp_genes =tpm %>% 
group_by(gene_id, gene_name, condition) %>% 
summarise(mean_tpm = mean(tpm)) %>% 
ungroup() %>% 
group_by(gene_id, gene_name) %>% 
summarise(max_mean_tpm = max(mean_tpm)) %>% 
ungroup() %>% 
filter(max_mean_tpm >= 1)

exp_genes %>% distinct(gene_name, gene_id)
# 34,026 genes left

# 
geneCount_exp = semi_join(
        geneCount,
        exp_genes
)

readCount_summary = left_join(
        read_delim(file.path(rnaseq_dir, "STAR","all.flagstat.txt"), delim="\t"),
        readCount_summary
)

write.csv(
        readCount_summary,
        file=file.path(rnaseq_dir,"HTseq/readCount_summary.csv"),
        row.names=F
)


dir.create(file.path(rnaseq_dir,"DE"), showWarnings=F)
save(geneCount, file=file.path(rnaseq_dir,"DE/geneCount.Rdata"))
save(tpm, file=file.path(rnaseq_dir,"DE/tpm.Rdata"))
save(meanTPM, file=file.path(rnaseq_dir,"DE/meanTPM.Rdata") )
save(exp_genes, file=file.path(rnaseq_dir,"DE/exp_genes.Rdata") )

########################### Jurkat ##############################
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/Jurkat2"

geneCount=read_delim(file.path(rnaseq_dir,"HTseq","htseq_counts.txt"), delim="\t")
geneCount=geneCount %>% filter(!grepl("__",gene))
colnames(geneCount)=gsub("_HTSeq","",colnames(geneCount))

# filter non-gencodeV19 genes
geneCount = geneCount %>% filter(grepl("ENSG", gene))

readCount_summary = tibble(
        sample=colnames(geneCount %>% select(-gene)),
        gene_align_count=apply(geneCount %>% select(-gene), 2, sum)
)

# read chrM and rRNA 
gene_id_rm = bind_rows(
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.rRNA_gene_id.txt", 
                delim="\t", col_names=c("gene_id","type")
        ),
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.chrM_gene_id.txt", 
                delim="\t", col_names=c("gene_id","type")
        )
) %>% distinct(gene_id)

# remove chrM and rRNA
geneCount = geneCount %>% filter(!gene %in% gene_id_rm$gene_id)

readCount_summary = left_join(
        readCount_summary,
        tibble(
                sample=colnames(geneCount %>% select(-gene)),
                gene_rm_chrM_rRNA_count=apply(geneCount %>% select(-gene), 2, sum)
        )
)

# remove genes using TPM < 1
gene_len=read_delim("/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.length.txt",delim="\t", 
col_names=c("gene","gene_name","gene_len"))

geneCount_melt = left_join(
        geneCount %>% 
        gather(key="sample",value="rawCount", -gene),
        gene_len
)
geneCount_melt = geneCount_melt %>% mutate(RPK=rawCount/gene_len*1000)

scale_factor=geneCount_melt %>% 
group_by(sample) %>% 
summarise(total_RPK=sum(RPK)) %>% 
ungroup() %>% 
mutate(scale_factor=total_RPK/1000000)

tpm=left_join(geneCount_melt, scale_factor) %>% 
mutate(tpm=RPK/scale_factor) %>% 
select(gene,gene_name,sample,tpm)

tpm = tpm %>% 
separate(sample, c("cell","time","ind"), sep="_", remove=F) %>% 
unite("condition",cell, time, sep="_")

meanTPM=tpm %>% 
group_by(gene, gene_name, condition) %>% 
summarise(tpm=mean(tpm)) %>% 
ungroup()

exp_genes =tpm %>% 
group_by(gene, gene_name, condition) %>% 
summarise(mean_tpm = mean(tpm)) %>% 
ungroup() %>% 
group_by(gene, gene_name) %>% 
summarise(max_mean_tpm = max(mean_tpm)) %>% 
ungroup() %>% 
filter(max_mean_tpm >= 1)

exp_genes %>% distinct(gene_name, gene)
# 26,326 genes left

# 
# geneCount_exp = semi_join(
#         geneCount,
#         exp_genes
# )

readCount_summary = left_join(
        read_delim(file.path(rnaseq_dir, "STAR","all.flagstat.txt"), delim="\t"),
        readCount_summary
)

write.csv(
        readCount_summary,
        file=file.path(rnaseq_dir,"HTseq/readCount_summary.csv"),
        row.names=F
)

dir.create(file.path(rnaseq_dir,"DE"), showWarnings=F)
save(geneCount, file=file.path(rnaseq_dir,"DE/geneCount.Rdata"))
save(tpm, file=file.path(rnaseq_dir,"DE/tpm.Rdata"))
save(meanTPM, file=file.path(rnaseq_dir,"DE/meanTPM.Rdata") )
save(exp_genes, file=file.path(rnaseq_dir,"DE/exp_genes.Rdata") )
