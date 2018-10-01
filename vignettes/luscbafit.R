library(iClusterPlus)
library(gplots)
library(lattice)
date()

## Load genomic datasets
gene = readRDS("~/iClusterPlus/data/gene.RDS")
seg.cn = readRDS("~/iClusterPlus/data/seg.cn.RDS")
meth450 = readRDS("~/iClusterPlus/data/meth450.RDS")
sample = readRDS("~/iClusterPlus/data/sample.RDS")
sample <- sample[!duplicated(sample$patient_id), ]
sample <- sample[order(sample$sample_id), ]

# get clinical data -----
lusc_mccv <- readRDS("~/iClusterPlusUtils/data/mc_lusc_clinical.RDS")
vfold <- lusc_mccv$splits[[1]] #choose splits
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
sample_id <- list( meth450$sample_id, rownames(seg.cn), sample$sample_id, gene$sample_id)
sample_id <- Reduce(intersect, sample_id)

## Filter
filter_sample_id <- function(s, m){
  m[rownames(m) %in% s, ]
}
## Filter
gene$matrix <- filter_sample_id(s = sample_id, m = gene$matrix)
meth450$matrix <- filter_sample_id(s = sample_id, m = meth450$matrix)
seg.cn <- filter_sample_id(s = sample_id, m = seg.cn)
sample <- sample[sample$sample_id %in% sample_id, ]

require(assertthat)
assert_that(all(rownames(meth450$matrix) == rownames(seg.cn)))
assert_that(all(rownames(gene$matrix) == rownames(seg.cn)))
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
                                              type = c("gaussian", "gaussian", "gaussian"),
                                              K= 1:8,
                                              n.burnin=10000,
                                              n.draw=12000,
                                              prior.gamma = rep(0.1, 6),
                                              sdev = 0.5,
                                              beta.var.scale = 1,
                                              thin = 1,
                                              pp.cutoff = 0.5)

saveRDS(tuned_fit, "~/iClusterPlus/lusc_ba/lusc.baclust_3dt.RDS")

