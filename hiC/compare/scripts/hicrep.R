# R/4.0.2
# remotes::install_github("aidenlab/straw/R")
# install.packages("/home/suc1/R/x86_64-pc-linux-gnu-library/4.0/hicrep_1.12.2.tar.gz", repo=NULL, type="source")
library(dplyr)
library(tidyr)
library(ggplot2)
library(hicrep)
library(gplots)
library(RColorBrewer)

args = commandArgs(trailingOnly=TRUE)
resol = as.integer(args[1])
out_prefix = args[2]
files = args[3:length(args)]

# resol=1000000
# out_prefix="/mnt/isilon/sfgi/suc1/analyses/grant/hiC/compare/pancreaticCells/plots/scc_heatmap.1M"
# files=dir(
#         "/mnt/isilon/sfgi/suc1/analyses/grant/hiC/hicup",
#         pattern="\\.1M\\.cool",
#         recursive=T,full.names=T
# )

mat_list=lapply(
        files,
        cool2matrix
)

# Resolution ----- h
# 10kb ----- 20
# 25kb ----- 10
# 40kb ----- 5
# 100kb ----- 3
# 500kb ----- 1 or 2
# 1Mb ----- 0 or 1

if(resol > 1000000){
        h=0
}else if(resol<=1000000 & resol>500000){
        h=1
}else if(resol<=500000 & resol>100000){
        h=2
}else if(resol<=100000 & resol>50000){
        h=3
}else if(resol<=50000 & resol>40000){
        h=4
}else if(resol<=40000 & resol>25000){
        h=5
}else if(resol<=25000 & resol>10000){
        h=10
}else if(resol<=10000){
        h=20
}


scc=matrix(nrow=length(mat_list),ncol=length(mat_list))

for(i in 1:(length(mat_list)-1)){
        total=i+1
        for(j in total:length(mat_list)){
                scc.out = get.scc(
                        mat_list[[i]], mat_list[[j]], 
                        resol = resol, h = h, 
                        lbr = 0, ubr = resol*100
                )
                scc[i,j]=scc.out$scc[1,1]
        }
}

for(i in 1:length(mat_list)){
        for(j in 1:length(mat_list)){
                if(i==j){
                        scc[i,j]=1
                }
                else if(j < i){
                        scc[i,j]=scc[j,i]
                }
        }
}

rownames(scc)= sapply(basename(files), function(x){strsplit(x,"\\.")[[1]][1]})
colnames(scc)= sapply(basename(files), function(x){strsplit(x,"\\.")[[1]][1]})

# plot heatmap
pdf(paste0(out_prefix,".pdf"))
p <- heatmap.2(
        scc,
        distfun = function(x) as.dist((1-x)/2),
        hclust=function(x) hclust(x,method="complete"),
        Rowv = TRUE, Colv= T,
        col=brewer.pal(9,"Blues"),
        symm = T,
        margins = c(12,12),
        trace="none",
        density.info="none",
        key.title = NA,
        keysize=1,
        key.xlab="SCC"
)
dev.off()

# write
scc_df = scc %>% as.data.frame() %>% 
tibble::rownames_to_column("sample1") %>% 
gather(key="sample2",value="scc",-sample1) %>% 
filter(sample1 > sample2)

write.csv(
        scc_df,
        file=paste0(out_prefix,".csv"), 
        row.names=F
)


