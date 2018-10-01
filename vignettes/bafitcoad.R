library(iClusterPlus)
library(gplots)
library(lattice)
date()

## Load genomic datasets
gene = readRDS("~/iClusterPlus/coad_results/gene_coad.RDS")
seg.cn = readRDS("~/iClusterPlus/coad_results/seg.cn_coad.RDS")
meth450 = readRDS("~/iClusterPlus/coad_results/meth450.RDS")
sample = readRDS("~/iClusterPlus/coad_results/sample_coad.RDS")
sample <- sample[!duplicated(sample$patient_id), ]
sample <- sample[order(sample$sample_id), ]
prot = readRDS("~/iClusterPlus/coad_results/prot.RDS")

# get clinical data -----
coad_mccv <- readRDS("~/iClusterPlusUtils/data/mc_coad_clinical.RDS")
vfold <- coad_mccv$splits[[1]] #choose splits
train_id <- vfold$in_id
clinical_data <- vfold$data

remove_is_na <- function(x){
  x[!is.na(x)]
}
is_in <- function(x, y){
  x %in% y
}

sample <- sample[is_in(sample$patient_id, clinical_data$patient_id), ]

## Get common sample id
sample_id <- list( meth450$sample_id, rownames(seg.cn), sample$sample_id, prot$sample_id, gene$sample_id)
sample_id <- Reduce(intersect, sample_id)

## Filter
filter_sample_id <- function(s, m){
  m[rownames(m) %in% s, ]
}
## Filter
gene$matrix <- filter_sample_id(s = sample_id, m = gene$matrix)
meth450$matrix <- filter_sample_id(s = sample_id, m = meth450$matrix)
seg.cn <- filter_sample_id(s = sample_id, m = seg.cn)
prot$matrix <- filter_sample_id(s = sample_id, m = prot$matrix )
sample <- sample[sample$sample_id %in% sample_id, ]

require(assertthat)
assert_that(all(rownames(meth450$matrix) == rownames(seg.cn)))
assert_that(all(rownames(gene$matrix) == rownames(seg.cn)))
assert_that(all(rownames(prot$matrix) == rownames(seg.cn)))
assert_that(all(sample$sample_id == rownames(seg.cn)))

# Filter clinical data
remove_is_na <- function(x){
  x[!is.na(x)]
}
clinical_data <- clinical_data[remove_is_na( match(sample$patient_id, clinical_data$patient_id, nomatch = NA_integer_)), ]


# get clinical data -----
remove_is_na <- function(x){
  x[!is.na(x)]
}
clinical_data <- clinical_data[remove_is_na( match(sample$patient_id, clinical_data$patient_id, nomatch = NA_integer_)), ]

tuned_fit <- iClusterPlus::tune.iClusterBayes(cpus = 8,
                                              dt1 = gene$matrix,
                                         dt2 = seg.cn,
                                         dt3 = meth450$matrix,
                                         dt4 = prot$matrix,
                                         type = c("gaussian","gaussian", "gaussian", "gaussian"),
                                         K= 1:8,
                                         n.burnin=10000,
                                         n.draw=12000,
                                         prior.gamma = rep(0.1, 6),
                                         sdev = 0.5,
                                         beta.var.scale = 1,
                                         thin = 1,
                                         pp.cutoff = 0.5)

saveRDS(tuned_fit, "~/iClusterPlus/coad_ba/coad.baclust_4dt.RDS")

