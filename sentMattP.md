# naive CD4 time-course

## Hi-C
what has been done
- pre-processing to matrix
`/mnt/isilon/sfgi/suc1/analyses/wells/hiC/hicup/naiveT_*`
- matrix merge, AB compartment, tads, loop calls
`/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_unstimulated_2reps`
`/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_8hr_2reps`
`/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_24hr_2reps`
- consensus loop calls
`/mnt/isilon/sfgi/suc1/analyses/wells/hiC/juicer/naiveT_stim_summary`
- global differential `/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/`
  - differential consensus loop `/mnt/isilon/sfgi/suc1/analyses/wells/hiC/compare/naiveT_stim/data/Rdata/sigDE_consensus_loop.Rdata`

## RNA-seq
`/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3`
- STAR
- HTseq
- DE
  - TPM `/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3/DE/tpm.Rdata`
  - pairwise differential `/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/PBMC_naiveT3/DE/DEG.Rdata`
  
## ATAC-seq
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT`
- encode 
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/CD4_*`
- peak: using pooledRep from each condition, filter by individual Rep peak (>= 2 reps), then merge condition to get consensus
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/peaks/mergeCondition_pooledRep/`
- annotation - promoter
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/annotation/promoter`
- bigwig
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/bigwig`
- quantitative differential
`/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/PBMC_naiveT/quantitative`


# NK, pDC
## ATAC
/mnt/isilon/sfgi/suc1/analyses/wells/atacSeq/CD_NK_cells
- encode

## RNAseq
`/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/CD_NK_cells`
- STAR
- HTseq
  - QC summary: `/mnt/isilon/sfgi/suc1/analyses/wells/rnaSeq/CD_NK_cells/HTseq/readCount_summary.csv`
  
