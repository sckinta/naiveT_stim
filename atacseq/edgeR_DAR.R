#R/4.0.2
if(!require("vsn")){
        BiocManager::install("vsn")
}
library(csaw)
library(edgeR)
library(statmod)
library(GenomicRanges)
library(RColorBrewer)
library(gplots)
library(tidyverse)
library(gridExtra)

dir="/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/quantitative"

load(file.path(dir,'all.filtered_newName.regCounts.Rdata'))
# load('param.RData')

############################ design_df ###################
design_df = design_df %>% 
select(sample, condition) %>% 
separate(sample, c("cell","time","rep"), sep="_", remove=F) %>% 
select(-cell,-time)

############################# PRE-FILTER ################################
colnames(y$counts) = design_df$sample
rownames(y$samples) = design_df$sample
y_cpm <- cpm(y) %>% as_tibble()
colnames(y_cpm) <- design_df$sample
y_cpm_transform <- y_cpm %>% 
mutate(id=rowRanges(reg.counts)$id) %>% 
gather(sample,cpm,-id)

dim(y$counts) # 181093      9

# top condition filter
total_sample=9
rep_sample=3
top_perc=1-(rep_sample-1)/total_sample
test <- y_cpm_transform %>% group_by(id) %>% summarise(q3=quantile(cpm,top_perc))
quantile(test$q3, seq(0.05,0.9,by=0.05))
# 5%       10%       15%       20%       25%       30%       35%       40% 
# 0.3354825 0.3928473 0.4353152 0.4785882 0.5198527 0.5638995 0.6131111 0.6694131 
# 45%       50%       55%       60%       65%       70%       75%       80% 
# 0.7350679 0.8175248 0.9274007 1.0712327 1.2753126 1.5713603 2.0038114 2.6957948 
# 85%       90% 
# 4.0015045 6.8414674 

sample_size = range(apply(y$counts, 2, sum)) # 5,119,691 16,616,337

cutoff=1
cutoff*sample_size/1000000 # means for each peak, we must have at least 5.119691 16.616337reads to support per peak
test %>% filter(q3>=cutoff) # 76,539 out of 181093 retained

y_cpm_transform_filtered <- semi_join(
        y_cpm_transform,
        test %>% filter(q3>=cutoff)
)
keep <- bind_rows(
        semi_join(y_cpm_transform ,y_cpm_transform_filtered ) %>% select(id) %>% mutate(keep=T),
        anti_join(y_cpm_transform ,y_cpm_transform_filtered ) %>% select(id) %>% mutate(keep=F)
) %>% separate(id,c("prefix",'id')) %>% select(id,keep)

keep <- keep %>% distinct() %>% arrange(as.numeric(id))
y_filtered <- y[unlist(keep %>% select(keep)),, keep.lib.sizes=FALSE]
colnames(y_filtered$counts) <- design_df$sample
rownames(y_filtered$samples) <- design_df$sample

sample_size_df = tibble(
        sample=design_df$sample,
        before=apply(y$counts, 2, sum),
        after=apply(y_filtered$counts, 2, sum)
) %>% mutate(retain_ratio=after/before)

range(sample_size_df$after) # 4587731 15199286
sample_size_df = sample_size_df %>% 
bind_rows(
        tibble(
                sample="peakN",
                before=nrow(y$counts),
                after=nrow(y_filtered$counts),
                retain_ratio=nrow(y_filtered$counts)/nrow(y$counts)
        )
)
write.csv(sample_size_df, file=file.path(dir,"sample_size_filter.csv"), row.names=F)


#############################  CPM #######################
y_cpm_filtered = cpm(y_filtered)
y_cpmVSN_filtered = normalizeVSN(y_cpm_filtered)
y_lcpm_filtered = cpm(y_filtered,log=T)

y_expAll_transform_filtered = data.frame(
        id=rowRanges(reg.counts)$id[keep$keep],
        id_num=rownames(y_cpm_filtered),
        y_cpm_filtered
) %>% as_tibble() %>% 
gather(key="sample", value="cpm", -id, -id_num) %>% 
left_join(
        data.frame(
                id=rowRanges(reg.counts)$id[keep$keep],
                id_num=rownames(y_cpmVSN_filtered),
                y_cpmVSN_filtered
        ) %>% as_tibble() %>% 
        gather(key="sample", value="cpmVSN", -id, -id_num)
) %>% 
left_join(
        data.frame(
                id=rowRanges(reg.counts)$id[keep$keep],
                id_num=rownames(y_lcpm_filtered),
                y_lcpm_filtered
        ) %>% as_tibble() %>% 
        gather(key="sample", value="lcpm", -id, -id_num)
)

precalled_peak_filtered = precalled_peak[keep$keep,]

save(y_expAll_transform_filtered, precalled_peak_filtered, y_filtered, design_df, keep, file=file.path(dir,"y_filtered.Rdata"))

############################### QC ###################################
dir.create(file.path(dir,"plots"), showWarnings=F)

### visualization of variation stablization
library(plyr)
sd2mean <- ddply(y_expAll_transform_filtered, .(id), summarize,
cpm_mean=mean(cpm),
cpm_sd=sd(cpm),
lcpm_mean=mean(lcpm),
lcpm_sd=sd(lcpm),
cpmVSN_mean=mean(cpmVSN),
cpmVSN_sd=sd(cpmVSN)
)
detach(package:plyr)
sd2mean <- sd2mean %>% as_tibble() %>% 
gather(type,value, -id) %>% 
separate(type, c("valueType","xy"), sep="_")
sd2mean <- sd2mean %>% spread(xy, value)

# plot sdVSmean
p1 <- ggplot(sd2mean,aes(mean,sd)) +
facet_wrap(~valueType, scale="free") +
geom_point() +
geom_smooth(color='red')
ggsave(file.path(dir,"plots/y_plotSDmean.png"), plot=p1, width = 21, height = 7)


### PCA
myPCAPlot <- function(df, type, design_df, shape_col, color_col, text_col=NULL) {
        myvars <- apply(as.matrix(df %>% select(id, sample, (!!type)) %>% spread(sample, (!!type)) %>% select(-id)), 1, var)
        d <- df %>% select(id, sample, (!!type)) %>% spread(sample, (!!type)) %>% mutate(sd=myvars) %>% filter(sd!=0) %>% select(-id, -sd)
        pca_m <- prcomp(t(d),
                         center = TRUE,
                         scale. = TRUE)
        pca_df <- data.frame(pca_m$x, sample=row.names(pca_m$x), stringsAsFactors = F) %>% tbl_df
        pca_df <- left_join(pca_df, design_df)
        eigs <- pca_m$sdev^2
        pca_summary <- cbind(
          SD = sqrt(eigs),
          Proportion = eigs/sum(eigs),
          Cumulative = cumsum(eigs)/sum(eigs)) %>% tbl_df
        pca_summary <- pca_summary %>%
          mutate(component=row_number())
        PC1_lab=unlist(pca_summary %>% 
                filter(component %in% c(1,2,3)) %>% 
                mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                filter(component==1) %>% 
                select(name),use.names=F
        )
        PC2_lab=unlist(pca_summary %>% 
                        filter(component %in% c(1,2,3)) %>% 
                        mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                        filter(component==2) %>% 
                        select(name),use.names=F
                )
        PC3_lab=unlist(pca_summary %>% 
                        filter(component %in% c(1,2,3)) %>% 
                        mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                        filter(component==3) %>% 
                        select(name),use.names=F
                )
                if (!is.null(shape_col)){
                        p1 <- ggplot(pca_df, aes(PC1,PC2)) +
                          geom_point(aes_string(col=color_col, shape=shape_col), size=5) +
                          xlab(PC1_lab) + ylab(PC2_lab) +
                          theme(panel.grid.minor=element_blank(), 
                          panel.grid.major=element_blank(),
                          panel.background=element_blank(),
                    panel.border=element_rect(colour = "black",fill=NA), 
                    axis.line=element_line())
                        #   geom_text(aes(label=reps), check_overlap = TRUE, size = 2)
                        p2 <- ggplot(pca_df, aes(PC1,PC3)) +
                          geom_point(aes_string(col=color_col, shape=shape_col), size=5) +
                          xlab(PC1_lab) + ylab(PC3_lab) +
                          theme(panel.grid.minor=element_blank(), 
                          panel.grid.major=element_blank(),
                          panel.background=element_blank(),
                    panel.border=element_rect(colour = "black",fill=NA), 
                    axis.line=element_line())
                        p3 <- ggplot(pca_df, aes(PC2,PC3)) +
                            geom_point(aes_string(col=color_col, shape=shape_col), size=5) +
                            xlab(PC2_lab) + ylab(PC3_lab) +
                            theme(panel.grid.minor=element_blank(), 
                          panel.grid.major=element_blank(),
                          panel.background=element_blank(),
                    panel.border=element_rect(colour = "black",fill=NA), 
                    axis.line=element_line())
            }else{
                    p1 <- ggplot(pca_df, aes(PC1,PC2)) +
                      geom_point(aes_string(col=color_col), size=5) +
                      xlab(PC1_lab) + ylab(PC2_lab) +
                      theme(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_blank(),
                      panel.background=element_blank(),
                panel.border=element_rect(colour = "black",fill=NA), 
                axis.line=element_line())
                    #   geom_text(aes(label=reps), check_overlap = TRUE, size = 2)
                    p2 <- ggplot(pca_df, aes(PC1,PC3)) +
                      geom_point(aes_string(col=color_col), size=5) +
                      xlab(PC1_lab) + ylab(PC3_lab) +
                      theme(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_blank(),
                      panel.background=element_blank(),
                panel.border=element_rect(colour = "black",fill=NA), 
                axis.line=element_line())
                    p3 <- ggplot(pca_df, aes(PC2,PC3)) +
                        geom_point(aes_string(col=color_col), size=5) +
                        xlab(PC2_lab) + ylab(PC3_lab) +
                        theme(panel.grid.minor=element_blank(), 
                      panel.grid.major=element_blank(),
                      panel.background=element_blank(),
                panel.border=element_rect(colour = "black",fill=NA), 
                axis.line=element_line())
            }
                
                if (!is.null(text_col)){
                    p1 <- p1 +
                      geom_text(aes_string(label=text_col))
                    p2 <- p2 +
                      geom_text(aes_string(label=text_col))
                    p3 <- p3 +
                      geom_text(aes_string(label=text_col))
            }
                
                p4 <- ggplot(pca_summary) +
                  geom_histogram(aes(component,Cumulative),
                                 stat="identity", fill="yellow", color="yellow", alpha=0.5) +
                  geom_hline(yintercept = 0.8, color="red") +
                  ylab('Cumulative proportion of variance') +
                  theme(panel.grid.minor=element_blank(), 
                  panel.grid.major=element_blank(),
                  panel.background=element_blank(),
            panel.border=element_rect(colour = "black",fill=NA), 
            axis.line=element_line())
                # grid.arrange(p1,p2,ncol=1)
                list(p1,p2,p3,p4,pca_df,pca_summary,PC1_lab,PC2_lab,PC3_lab)
}

p <- myPCAPlot(
        y_expAll_transform_filtered, 
        "cpmVSN",
        design_df=design_df,
        shape_col="rep", # design_df
        color_col="condition", # design_df
        text_col=NULL
)
pdf(file.path(dir,"plots/PCA_3_allOCRs.pdf"),width=21, height=7)
do.call("grid.arrange", c(p[1:3], ncol=3))
dev.off()

### heatmap
myHeatmap <- function (df) {
  sample_cor <- as.matrix(cor(df))
  heatmap.2(sample_cor,
            distfun = function(x) as.dist((1-x)/2),
            hclust=function(x) hclust(x,method="average"),
            Rowv = TRUE, Colv= T,
            col=brewer.pal(9,"Blues"),
            symm = T,
            margins = c(12,12),
            trace="none",
            density.info="none",
            key.title = NA,
            keysize=1,
            key.xlab="correlation")
}

mysds <- apply(
        as.matrix(
                y_expAll_transform_filtered %>% 
                select(id, sample, cpmVSN) %>% 
                spread(sample, cpmVSN) %>% 
                select(-id)
        )
        , 1, sd
)
d <- y_expAll_transform_filtered %>% 
select(id, sample, cpmVSN) %>% 
spread(sample, cpmVSN) %>% 
mutate(sd=mysds) %>% 
filter(sd!=0) %>% 
select(-id, -sd)

pdf(file.path(dir,'plots/y_cpmVSN_correlation_heatmap.pdf'))
myHeatmap(d)
dev.off()

rm(list=ls())

############################ pairwise DEA ######################################
load(file.path(dir,"y_filtered.Rdata"))

dir.create(file.path(dir,"DE"))
dir.create(file.path(dir,"DE/tables"))

## design matrix
group = design_df$condition
group = factor(group, levels=c("CD4_unstim", "CD4_8hr", "CD4_24hr"))
rep = design_df$rep
design <- model.matrix(~rep+group)
colnames(design) <- gsub("group|rep","",colnames(design))

## estimateDisp
y_filtered <- estimateDisp(y_filtered,design)

## fit
fit <- glmQLFit(y_filtered, design, robust=T)

## contrasts
my.contrasts <- makeContrasts(
        CD4_8hr_vs_unstim = CD4_8hr,
        CD4_24hr_vs_unstim = CD4_24hr,
        CD4_24hr_vs_8hr = CD4_24hr - CD4_8hr,
        levels=design
)

## run pairwise comparison glmQLFTest in loop
writeTest <- function(name,dir){
        results = glmQLFTest(fit,contrast=my.contrasts[,name])
        FDR = p.adjust(results$table$PValue, 'fdr')
        
        df = bind_cols(
                precalled_peak_filtered %>% 
                as.data.frame() %>% tbl_df %>% 
                select(id, chr=seqnames, start, end),
                data.frame(
                        results$table, 
                        FDR=FDR
                ) %>% tbl_df
        ) %>% arrange(FDR)
        filename = paste0(name,".DE_edgeR.txt")
        write.table(df, file=file.path(dir,filename), sep="\t", row.name=F, quote=F)
        df %>% mutate(comp=name)
}

DE_list = lapply(colnames(my.contrasts),writeTest,dir=file.path(dir,"DE/tables"))
DE_df = do.call("bind_rows",DE_list)

## sigDE
sigDE_df = DE_df %>% 
filter(FDR < 0.05, abs(logFC)>1) %>% 
mutate(direction=ifelse(logFC > 1,"up","down"))

sigDE_df %>% distinct(id)
# 35,926 OCRs were sigDE among all pairwise

writeSigDE <- function(name, dir){
        df=sigDE_df %>% filter(comp==name) %>% select(-comp)
        filename=paste0(name,".sigDE_edgeR.txt")
        write.table(df, file=file.path(dir, filename), sep="\t", row.names=F, quote=F)
}

lapply(colnames(my.contrasts),writeSigDE, dir=file.path(dir,"DE/tables"))

### summarize DAR
out = DE_df %>% 
filter(FDR < 0.05, abs(logFC)>1) %>% 
mutate(direction=ifelse(logFC > 1,"up","down")) %>% 
group_by(direction,comp) %>% 
summarise(OCR_n=n_distinct(id)) %>% 
ungroup()

save(y_filtered, DE_df, sigDE_df, out, my.contrasts, file=file.path(dir,"DE","DEA.Rdata"))


############################### annotate DAR to promoter ########################
prOCR2gene = read_delim("/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/annotation/promoter/prOCR2gene.txt", delim="\t")

# promoter OCRs
prOCR_sigDE = sigDE_df %>% 
distinct(id, direction, comp) %>% 
semi_join(
        prOCR2gene %>% 
        select(id=ocr_id, gene_name, gene_id, gene_type) %>% 
        distinct()
) %>% left_join(
        prOCR2gene %>% 
        select(id=ocr_id, gene_name, gene_id, gene_type) %>% 
        distinct()
)

prOCR_sigDE %>% 
filter(gene_type != "pseudogene") %>% 
summarise(
        OCR_n=n_distinct(id),
        gene_n=n_distinct(gene_id)
)# 8062   DAR corresponding to  7442 nonpseudogene promoter


out2 =  prOCR_sigDE %>% 
filter(gene_type != "pseudogene") %>% 
group_by(direction,comp) %>% 
summarise(
        prOCR_n = n_distinct(id),
        gene_n = n_distinct(gene_id)
) %>% ungroup()


out = bind_rows(
        left_join(out, out2) %>% 
        arrange(comp),
        sigDE_df %>% 
        group_by(comp) %>% 
        summarise(OCR_n = n_distinct(id)) %>% 
        left_join(
                prOCR_sigDE %>% 
                filter(gene_type != "pseudogene") %>% 
                group_by(comp) %>% 
                summarise(
                        prOCR_n = n_distinct(id),
                        gene_n = n_distinct(gene_id)
                )
        ) %>% mutate(direction="total")
) %>% arrange(comp)

write.csv(out, file=file.path(dir,"DE/tables","sigDE_summary.csv"), row.names=F)

######################## FPKM #########################
load(file.path(dir,'all.filtered_newName.regCounts.Rdata'))
design_df = design_df %>% 
select(sample, condition)

# fpkm
gr = rowRanges(reg.counts)
counts=y$counts
colnames(counts) = design_df %>% pull(sample)

df = data.frame(id_num=rownames(counts),
id=gr$id,
chr=seqnames(gr),
start=start(gr),
end=end(gr),
counts,
peak_len = end(gr) - (start(gr)-1)) %>% tbl_df

df = df %>% gather(key="sample",value="count",-id_num,-id,-chr,-start,-end,-peak_len)

lib.size=data.frame(sample=colnames(counts),lib.size=colSums(y$counts)) %>% tbl_df
df=left_join(df,lib.size)

df = df %>% mutate(fpkm=count/(lib.size*peak_len)*10^9) %>% 
mutate(log2fpkm=log2((count+1)/(lib.size*peak_len)*10^9))
fpkm = df %>% select(-id_num, -lib.size, -peak_len)
fpkm = fpkm %>% 
left_join(design_df %>% select(sample, condition))
save(fpkm , file=file.path(dir,"fpkm.Rdata"))

### meanFpkm
meanFpkm = left_join(
        fpkm,
        design_df
) %>% 
group_by(id, chr, start, end, condition) %>% 
summarise(fpkm=mean(fpkm), log2fpkm=mean(log2fpkm)) %>% 
ungroup()

save(meanFpkm , file=file.path(dir,"meanFpkm.Rdata"))
