# R/4.0.2
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

dir="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells"

######################### pileup across samples ########################
files=dir(file.path(dir,"data/40K/eigen"), pattern="cis\\.vecs\\.tsv", full.names=T) # 464424 lines

df = lapply(
        files,
        function(file){
                sample=strsplit(basename(file),"\\.")[[1]][1]
                read_delim(file, delim="\t") %>% 
                select(chrom, start, end, GC, E1, E2) %>% 
                mutate(sample=sample)
        }
)

df = do.call("bind_rows", df)

eigen_E1 = df %>% filter(!is.na(E1)) %>% 
select(-E2) %>% 
spread(key=sample, value=E1)

eigen_E2 = df %>% filter(!is.na(E2)) %>% 
select(-E1) %>% 
spread(key=sample, value=E2)

save(eigen_E1, file=file.path(dir,"data/40K/eigen/eigen_E1.Rdata"))
save(eigen_E2, file=file.path(dir,"data/40K/eigen/eigen_E2.Rdata"))


######################## call compartments ###############################
files=dir(file.path(dir,"data/40K/eigen"), pattern="cis\\.vecs\\.tsv",full.names=T) # 464424 lines

df = lapply(
        files,
        function(file){
                sample=strsplit(basename(file),"\\.")[[1]][1]
                read_delim(file, delim="\t") %>% 
                select(chrom, start, end, GC, E1, E2) %>% 
                mutate(sample=sample)
        }
)

tmp = df[[1]] %>% 
filter(!is.na(E1)) %>% 
mutate(E1_sign=sign(E1))

tmp_list=split(tmp, tmp$chrom)
chr_tmp = tmp_list[["chr20"]]

chr_tmp = chr_tmp %>% 
mutate(lag_sign=lag(chr_tmp$E1_sign)) %>%
mutate(break_here=ifelse(E1_sign*lag_sign==-1, 1,0))

which(chr_tmp$break_here==1)
