library(tidyverse)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(gplots)
library(ggrepel)

rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3"
setwd(file.path(rnaseq_dir,"DE"))

################### DESIGN ######################
load(file.path(rnaseq_dir,"DE/geneCount.Rdata"))
exp_design = tibble(
        sample=colnames(geneCount %>% select(-gene_id, -gene_name, -gene_type, -biotype))
) %>% 
separate(sample, c("cell","time","rep"), sep="_", remove=F) %>% 
mutate(condition=glue::glue("{cell}_{time}"))

group=factor(exp_design$condition)

######### DEGList FILTER LOW EXPRESSION GENES ############
y=DGEList(
        counts=geneCount %>% 
        select(-gene_name, -gene_type, -biotype) %>% 
        as.data.frame() %>% 
        column_to_rownames("gene_id"),
        group=group
)

y_cpm = cpm(y)
nrow(y$samples) # 9
dim(y$counts)  # 41715     9
range(y$samples$lib.size)  # 16212497 58523513 (16M to 58M)
top_sample_perc=1-(3-1)/9 ### assume top 3 out of 9 are from the same condition
quantile_cpm=apply(y_cpm,1,function(x){quantile(x,top_sample_perc)})
quantile_cpm_quantile=quantile(
        quantile_cpm,
        c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8)
)
quantile_cpm_quantile

# 10%        15%        20%        25%        30%        35%        40% 
# 0.1244102  0.2043951  0.3015504  0.4059736  0.5245369  0.6735446  0.8686683 
#  50%        60%        70%        80% 
# 1.5773464  2.9611994  6.3875938 20.9010554

cpm_cutoff=0.65  #### filter out lowest 60-70% genes
range(y$samples$lib.size)/10^6*cpm_cutoff  # assume you need as little as ~ min_count reads in smallest samples to prove that gene express
# 10.53812 38.04028
sum(quantile_cpm > cpm_cutoff) # 27404 genes remain for differential analysis

# filter now
keep=(quantile_cpm > cpm_cutoff)
y <- y[keep, , keep.lib.sizes=FALSE]
range(y$samples$lib.size)  # 16145481 58503412 (16M to 58.5M)
dim(y$counts) # 27404     9


####################### y_expAll_melt ###################
y <- calcNormFactors(y,method="TMM")

y_cpm=data.frame(gene=rownames(cpm(y)),cpm(y)) %>% as_tibble() #use normalized lib.size read_count/(norm.factor*lib.size)

y_cpm_log=data.frame(gene=rownames(cpm(y,log=T)),cpm(y,log=T)) %>% as_tibble() 

y_cpmVSN=normalizeVSN(y)
y_cpmVSN=data.frame(gene=rownames(y_cpmVSN),y_cpmVSN) %>% as_tibble() 

y_expAll_melt=bind_rows(
        y_cpm %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpm"),
        y_cpm_log %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpm_log"),
        y_cpmVSN %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpmVSN")
) %>% arrange(gene)
y_expAll_melt=y_expAll_melt %>% spread(exp_type,exp)

save(y, y_expAll_melt, exp_design, keep,file=file.path(rnaseq_dir,"DE","y_data.RData"))

####################### QC plot #############################
### MSDplot
library(plyr)
sd2mean <- ddply(y_expAll_melt, .(gene), summarize,
cpm_mean=mean(cpm),
cpm_sd=sd(cpm),
lcpm_mean=mean(cpm_log),
lcpm_sd=sd(cpm_log),
cpmVSN_mean=mean(cpmVSN),
cpmVSN_sd=sd(cpmVSN)
)
detach(package:plyr)

sd2mean <- sd2mean %>% tbl_df %>% gather(type,value, -gene) %>% separate(type, c("valueType","xy"), sep="_")
sd2mean <- sd2mean %>% spread(xy, value)

p1 <- ggplot(sd2mean,aes(mean,sd)) +
facet_wrap(~valueType, scale="free") +
geom_point() +
geom_smooth(color='red')
dir.create(file.path(rnaseq_dir, "DE/plots"), showWarnings = FALSE)
ggsave(file.path(rnaseq_dir,"DE/plots/y_plotSDmean.png"), plot=p1, width = 21, height = 7)

### PCA plot
myPCAPlot <- function(df, exp_design, type, shape_col, color_col, text_col=NULL) {
        myvars <- apply(as.matrix(df %>% select(gene, sample, (!!type)) %>% spread(sample, (!!type)) %>% select(-gene)), 1, var)
        d <- df %>% select(gene, sample, (!!type)) %>% spread(sample, (!!type)) %>% mutate(sd=myvars) %>% filter(sd!=0) %>% select(-gene, -sd)
        pca_m <- prcomp(t(d),
                         center = TRUE,
                         scale. = TRUE)
        pca_df <- data.frame(pca_m$x, sample=row.names(pca_m$x), stringsAsFactors = F) %>% tbl_df
        pca_df <- left_join(pca_df, exp_design)
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
                      geom_text_repel(aes_string(label=text_col))
                    p2 <- p2 +
                      geom_text_repel(aes_string(label=text_col))
                    p3 <- p3 +
                      geom_text_repel(aes_string(label=text_col))
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
        y_expAll_melt, 
        exp_design,
        "cpmVSN",
        shape_col="rep", # exp_design
        color_col="condition", # exp_design
        text_col=NULL
)
pdf(file.path(rnaseq_dir,"DE/plots/PCA_3.pdf"),width=21, height=7)
do.call("grid.arrange", c(p[1:3], ncol=3))
dev.off()

### pairwise heatmap
myHeatmap <- function (df) {
  sample_cor <- as.matrix(cor(df))
  heatmap.2(sample_cor,
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
            key.xlab="correlation")
}

mysds <- apply(as.matrix(y_expAll_melt %>% select(gene, sample, cpmVSN) %>% spread(sample, cpmVSN) %>% select(-gene)), 1, sd)
d <- y_expAll_melt %>% select(gene, sample, cpmVSN) %>% spread(sample, cpmVSN) %>% mutate(sd=mysds) %>% filter(sd!=0) %>% select(-gene, -sd)
pdf(file.path(rnaseq_dir,"DE/plots/corr_heatmap.pdf"))
myHeatmap(d)
dev.off()

#################### pairwise-DEG ##############################
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3"
load(file.path(rnaseq_dir,"DE/y_data.RData"))

## design matrix
group = factor(exp_design$condition, levels=c("CD4_unstim","CD4_8hr","CD4_24hr")
rep = exp_design$rep
design <- model.matrix(~rep+group)
colnames(design) <- gsub("group|rep","",colnames(design))

## estimateDisp
y <- estimateDisp(y,design)

## fit
fit <- glmQLFit(y, design, robust=T)

## contrasts
my.contrasts <- makeContrasts(
        CD4_24hr_vs_unstim=CD4_24hr,
        CD4_8hr_vs_unstim=CD4_8hr,
        CD4_24hr_vs_8hr=CD4_24hr-CD4_8hr,
        levels=design
)

## run pairwise comparison glmQLFTest in loop
gene_type = read_delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode.v19.gene_type.txt", delim="\t", comment="#")

writeTest <- function(name,dir){
        results = glmQLFTest(fit,contrast=my.contrasts[,name])
        FDR = p.adjust(results$table$PValue, 'fdr')
        df = data.frame(
                gene_id=rownames(results$table),
                results$table, FDR=FDR
        ) %>% tbl_df %>% arrange(FDR)
        df = left_join(
                df,
                gene_type %>% distinct(gene_id, gene_name)
        ) %>% select(gene_name, gene_id, logFC, logCPM, F, PValue, FDR)
        filename = paste0(name,".DE_edgeR.txt")
        write.table(df, file=file.path(dir,filename), sep="\t", row.name=F, quote=F)
        df %>% mutate(comp=name)
}

dir.create(file.path(rnaseq_dir,"DE/tables"),showWarnings = FALSE)
DE_list = lapply(colnames(my.contrasts),writeTest,dir=file.path(rnaseq_dir,"DE/tables"))
DE_df = do.call("bind_rows",DE_list)

## sigDE
out = DE_df %>% 
filter(FDR < 0.05, abs(logFC)>1) %>% 
mutate(direction=ifelse(logFC > 1,"up","down")) %>% 
group_by(direction,comp) %>% summarise(DEgene_num=n_distinct(gene_id)) %>% 
spread(key=direction,value=DEgene_num) %>% mutate(total=down+up)

write.csv(out, file=file.path(rnaseq_dir,"DE/tables","sigDE_summary.csv"), row.names=F)

sigDE_df = DE_df %>% 
filter(FDR < 0.05, abs(logFC)>1) %>% 
mutate(direction=ifelse(logFC > 1,"up","down"))

sigDE_df %>% distinct(gene_id)
# 9,815 genes were sigDE among all pairwise

writeSigDE <- function(name, dir){
        df=sigDE_df %>% filter(comp==name) %>% select(-comp)
        filename=paste0(name,".sigDE_edgeR.txt")
        write.table(df, file=file.path(dir, filename), sep="\t", row.names=F, quote=F)
}

lapply(colnames(my.contrasts),writeSigDE, dir=file.path(rnaseq_dir,"DE/tables"))

save(y, DE_df, sigDE_df, out, file=file.path(rnaseq_dir,"DE","DEG.Rdata"))



###################### downstream analysis ########################
rnaseq_dir="/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3"
load(file.path(rnaseq_dir,"DE","DEG.Rdata"))
load(file.path(rnaseq_dir,"DE","tpm.Rdata"))
load(file.path(rnaseq_dir,"DE","meanTPM.Rdata"))
load(file.path(rnaseq_dir,"DE","exp_genes.Rdata"))

DEG_log2tpm = semi_join(
        tpm,
        sigDE_df %>% 
        distinct(gene_id)
) %>% 
semi_join(exp_genes) %>% 
mutate(log2tpm=log2(tpm+1)) %>% 
select(gene_name, gene_id, sample, log2tpm) %>% 
spread(key=sample, value=log2tpm) # 9,525 genes with max_mean_tpm > 1

DEG_log2tpm

# prepare matrix to plot
DEG_log2tpm_m = DEG_log2tpm %>% 
select(gene_id, contains("unstim"), contains("8hr"), contains("24hr")) %>% 
as.data.frame() %>% 
column_to_rownames("gene_id") %>% 
as.matrix()

DEG_log2tpm_scale_m=t(apply(DEG_log2tpm_m, 1, scale))
colnames(DEG_log2tpm_scale_m) = colnames(DEG_log2tpm_m)

DEG_meanlog2tpm_m = apply(DEG_log2tpm_m,1,mean)

DEG_log10FDR_m = DE_df %>% 
semi_join(DEG_log2tpm) %>% 
mutate(log10FDR=-log10(FDR)) %>% 
select(comp, gene_id, log10FDR) %>% 
spread(key=comp, value=log10FDR) %>% 
select(gene_id, CD4_8hr_vs_unstim, CD4_24hr_vs_unstim, CD4_24hr_vs_8hr) %>% 
as.data.frame() %>% 
column_to_rownames("gene_id") %>% 
as.matrix()
DEG_log10FDR_m = DEG_log10FDR_m[names(DEG_meanlog2tpm_m),]

#### plot nosplit
library(ComplexHeatmap)
library(circlize)

p = Heatmap(
        DEG_log2tpm_scale_m,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"),space="RGB"), 
        name = "Expression (log2TPM) Z-score", 
        cluster_rows = T,
        clustering_distance_rows=function(x) as.dist(1-cor(t(x))),
        clustering_method_rows="average",
        show_row_names = F,
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(5, "cm")
) +
Heatmap(
        DEG_meanlog2tpm_m,
        name="Avg. Expression (log2TPM)",
        col = colorRamp2(c(0, 10), c("white", "blue"),space="RGB"),  # range(DEG_meanlog2tpm_m)
        width = unit(5, "mm")
) +
Heatmap(
        DEG_log10FDR_m,
        name="DEG -log10FDR",
        col = colorRamp2(c(0, 10), c("white", "orange"),space="RGB"),
        show_row_names = F,
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(15, "mm")
)

X11.options(colortype="pseudo.cube")
pdf(file.path(rnaseq_dir,"DE/plots","DEG_tpm_heatmap.pdf"), width=10)
draw(p)
dev.off()
X11.options(reset = T)

### hierachical split
# cutree(hclust(dd, method="average"), k = 5)
library(dendextend)
dd = as.dist(1-cor(t(DEG_log2tpm_scale_m)))
hc = hclust(dd, method="average")
pdf(file.path(rnaseq_dir,"DE/plots","tmp.pdf"))
plot(hc)
dev.off()
row_dend = as.dendrogram(hc)

# pdf(file.path(rnaseq_dir,"DE/plots","tmp.pdf"))
# row_dend %>% 
# set("labels", NULL) %>% 
# set("branches_k_color", value = 1:5, k=5) %>% 
# plot()
# dev.off()

row_dend = color_branches(row_dend, k = 3)
# table(cutree(row_dend, k = 3))


p = Heatmap(
        DEG_log2tpm_scale_m,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"),space="RGB"), 
        name = "Expression (TPM) Z-score", 
        cluster_rows = row_dend,
        row_gap=unit(1,"mm"),
        show_row_names = F,
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(5, "cm")
) +
Heatmap(
        DEG_meanlog2tpm_m,
        name="Avg. Expression (log2TPM)",
        col = colorRamp2(c(0, 10), c("white", "blue"),space="RGB"),  # range(DEG_meanlog2tpm_m)
        width = unit(5, "mm")
) +
Heatmap(
        DEG_log10FDR_m,
        name="DEG -log10FDR",
        col = colorRamp2(c(0, 10), c("white", "orange"),space="RGB"),
        show_row_names = F,
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(15, "mm")
)

X11.options(colortype="pseudo.cube")
pdf(file.path(rnaseq_dir,"DE/plots","DEG_tpm_heatmap_split.pdf"), width=10)
draw(p)
dev.off()
X11.options(reset = T)

### kmeans split
set.seed(123)
km = kmeans(as.dist(1-cor(t(DEG_log2tpm_scale_m))), centers = 5)
table(km$cluster)

p = Heatmap(
        DEG_log2tpm_scale_m,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"),space="RGB"), 
        name = "Expression Z-score", 
        cluster_rows = T,
        clustering_distance_rows=function(x) as.dist(1-cor(t(x))),
        clustering_method_rows="average",
        show_row_names = F,
        row_km = 5,
        row_title_gp = gpar(fill = c("red", "blue", "green", "orange", "purple"), fontsize=rep(10,5)),
        gap = unit(1, "mm"),
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(5, "cm")
) +
Heatmap(
        DEG_meanlog2tpm_m,
        name="Avg. Expression (log2TPM)",
        col = colorRamp2(c(0, 10), c("white", "blue"),space="RGB"),  # range(DEG_meanlog2tpm_m)
        width = unit(5, "mm")
) +
Heatmap(
        DEG_log10FDR_m,
        name="DEG -log10FDR",
        col = colorRamp2(c(0, 10), c("white", "orange"),space="RGB"),
        show_row_names = F,
        cluster_columns= F,
        show_column_dend= F,
        show_column_names = T, 
        width = unit(15, "mm")
)

X11.options(colortype="pseudo.cube")
pdf(file.path(rnaseq_dir,"DE/plots","DEG_tpm_heatmap_split.pdf"), width=10)
set.seed(123)
draw(p)
dev.off()
X11.options(reset = T)

