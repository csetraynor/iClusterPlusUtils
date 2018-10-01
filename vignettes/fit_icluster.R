library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
source("R/preproc.R")
date()

gene = readRDS("data/gene.RDS")
seg.cn = readRDS("data/seg.cn.RDS")
meth450 = readRDS("data/meth450.RDS")
sample = readRDS("data/sample.RDS")

## Get common sample id
sample_id <- list( gene$sample_id, rownames(seg.cn), sample$sample_id, meth450$sample_id)
sample_id <- Reduce(intersect, sample_id)

## Filter
gene$matrix <- filter_sample_id(s = sample_id, m = gene$matrix)
meth450$matrix <- filter_sample_id(s = sample_id, m = meth450$matrix)
seg.cn <- filter_sample_id(s = sample_id, m = seg.cn)


require(assertthat)
assert_that(all(rownames(meth450$matrix) == rownames(seg.cn)) )
assert_that(all(rownames(gene$matrix) == rownames(seg.cn)) )

out_folder <- "lusc_icluster_results/"

set.seed(321)
for(k in 6:20){
  cv.fit = tune.iClusterPlus(cpus = 28,
                             dt1 = gene$matrix,
                             dt2 = seg.cn,
                             dt3 = meth450$matrix,
                             type = c("gaussian","gaussian", "gaussian"),
                             n.lambda= NULL, K=k, maxiter=20)
  saveRDS(cv.fit, file = paste0(out_folder, "cv.fit_3.k",k,".RDS"))
}
# 
date()
