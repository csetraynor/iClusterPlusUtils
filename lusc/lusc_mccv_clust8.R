#mc-cross-validation.R
# Load libraries
suppressMessages(library(dplyr))
# For glmnet
suppressMessages(library(glmnet))
suppressMessages(library(glmnetUtils))
# For iCluster
library(iClusterPlus)
library(gplots)
library(lattice)
# For survival analysis
library(survival)
library(survcomp)
# For MC-crossvalidations
suppressMessages(library(rsample))
# for this package
#Read arguments listed on command line
args = (commandArgs(TRUE))
i = as.integer(args[1]) 
print(i)
set.seed(69837693)

## Load genomic datasets
gene = readRDS("~/iClusterPlus/data/gene.RDS")
seg.cn = readRDS("~/iClusterPlus/data/seg.cn.RDS")
meth450 = readRDS("~/iClusterPlus/data/meth450.RDS")
sample = readRDS("~/iClusterPlus/data/sample.RDS")
## Get common sample id
sample_id <- list(rownames(gene$matrix), rownames(seg.cn), sample$sample_id, 
                  rownames(meth450$matrix))
sample_id <- Reduce(intersect, sample_id)

## Filter
filter_sample_id <- function(s, m){
  m[rownames(m) %in% s, ]
}
gene$matrix <- filter_sample_id(s = sample_id, m = gene$matrix)
meth450$matrix <- filter_sample_id(s = sample_id, m = meth450$matrix)
seg.cn <- filter_sample_id(s = sample_id, m = seg.cn)
sample <- sample[sample$sample_id %in% sample_id, ]
sample <- sample[order(sample$sample_id), ]

require(assertthat)
assert_that(all(rownames(meth450$matrix) == rownames(seg.cn)))
assert_that(all(rownames(gene$matrix) == rownames(seg.cn)))
## Select sigfeatures

# get clinical data -----
lusc_mccv <- readRDS("~/iClusterPlusUtils/data/mc_lusc_clinical.RDS")
vfold <- lusc_mccv$splits[[i]] #choose splits
train_id <- vfold$in_id
date()
### --- Conduct integrative Clustering
set.seed(321)
model.tuned <- iClusterPlus::tune.iClusterBayes(cpus = 8,
                                                dt1 = gene$matrix[train_id, ],
                                                dt2 = seg.cn[train_id, ],
                                                dt3 = meth450$matrix[train_id, ],
                                                type = c("gaussian","gaussian", "gaussian"),
                                                K= 1:8,
                                                n.burnin=10000,
                                                n.draw=12000,
                                                prior.gamma = rep(0.1, 6),
                                                sdev = 0.5,
                                                beta.var.scale = 1,
                                                thin = 1,
                                                pp.cutoff = 0.5)
date()
## Model selection
nK = length(model.tuned$fit)

## get the optimal k+1 or g cluster model
m.BIC <- sapply(model.tuned$fit, function(m){
  m$BIC
})
k = which.min(m.BIC)
best.fit <- model.tuned$fit[[k]]

model.test <- iClusterPlus::iClusterBayes(dt1 = gene$matrix[-train_id, ],
                                          dt2 = seg.cn[-train_id, ],
                                          dt3 = meth450$matrix[-train_id, ],
                                          type = c("gaussian","gaussian", "gaussian"),
                                          K= k,
                                          n.burnin=1000,
                                          n.draw=1200,
                                          prior.gamma = rep(0.1, 6),
                                          sdev = 0.5,
                                          beta.var.scale = 1,
                                          thin = 1,
                                          pp.cutoff = 0.5)


cluster_predict_kmeans <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  max.col(-t(tmp))  # find index of min distance
}

centers <- best.fit$centers
x <- model.test$meanZ


## Include cluster variable in model
train <- rsample::analysis(vfold)
test <- rsample::assessment(vfold)

train$cluster <- as.factor( best.fit$clusters )
test$cluster <- as.factor( cluster_predict_kmeans(x, centers) )

# Free memory
rm(list = c("gene", "seg.cn", "meth450", "model.test", "model.tuned"))

### Standardise covariates ----------------------
preProcValues <- caret::preProcess(train[ "age" ], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, train[ "age" ] ) )
colnames(trainTransformed) <- "age"

#colnames(trainTransformed) = c("age_std")
train <- cbind(train[ ,-match(c("age"), colnames(train) )], trainTransformed) %>%
  dplyr::mutate(time = os_months,
                status = os_status)
rm(trainTransformed) # but we keep preProcValues
## transform test set with training vars
testTransformed <- as.data.frame( predict(preProcValues, test[ "age" ]) )
test <- cbind(test[ ,-match(c("age"), colnames(test) )], testTransformed) %>%
  dplyr::mutate(time = os_months,
                status = os_status)

## Cox model
clinical_vars <- c("stage", "smoke", "age")
surv_formula_clinical <- as.formula(paste0("Surv(time, status) ~ ", 
                                           paste(clinical_vars, collapse = "+")))
## Create model matrix --
x_clinical <- model.matrix(surv_formula_clinical, data = train)[, 
                                                                -1]  #drop the intercept
x_clinical_test <- model.matrix(surv_formula_clinical, data = test)[, 
                                                                    -1]  #drop the intercept
cox.train <- coxph(surv_formula_clinical, data = train)
# calculate risk score predictions test and train

cox_rp = as.matrix(x_clinical_test) %*% as.vector(unlist(coef(cox.train)))

cox_rp_train = as.matrix(x_clinical) %*% as.vector(unlist(coef(cox.train)))




# calculate metrics
get_survmetrics <- function(train, test = NULL, cox_rp,  cox_rp_train = NULL){
  
  if(is.null(test)){
    test <- train
    cox_rp_train <-  cox_rp
  }
  
  ## Prepare to calculate metrics
  ws <- rep(1, length(test$time))
  ws[test$time > max(test$time[test$status] ) ] <- 0.01
  
  obs.test.time <- sort ( unique( test$time[test$status] ) )
  #cindex
  Cindex.cox <- survcomp::concordance.index(x = cox_rp, surv.time = test$time, surv.event = test$status, method = "conservative",  na.rm=TRUE, weights = ws)$c.index
  
  #brier score
  bs.cox <- survcomp::sbrier.score2proba(data.tr = data.frame(time = train$time, event = train$status, score = cox_rp_train), data.ts = data.frame(time = test$time, event = test$status, score = cox_rp), method = c("cox"))
  
  #dindex
  D.cox <- survcomp::D.index(x = cox_rp, surv.time = test$time, surv.event = test$status, weights = ws)
  
  #AUC ROC
  AUC.cox <- survcomp::tdrocc(x = as.vector(cox_rp) , surv.time = test$time, surv.event = test$status,  time = max(obs.test.time) )
  
  list("Cindex.cox" = Cindex.cox,  "bs.cox" = bs.cox, "AUC.cox" = AUC.cox, "D.cox" = D.cox) 
}

metrics_list <- get_survmetrics(train = train, test = test, cox_rp = cox_rp, 
                                cox_rp_train = cox_rp_train)
cindex <- c(metrics_list[["Cindex.cox"]])
bs <- c(metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(metrics_list[["AUC.cox"]]$AUC)
D <- c(metrics_list[["D.cox"]]$d.index)
models <- c("Cox's-clinical")
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc)
## Genomic models
surv_formula <- as.formula(paste0("Surv(time, status) ~ cluster") )

## Create model matrix --
x <- model.matrix(surv_formula, data = train)[ ,-1]  #drop the intercept
x_test <- model.matrix(surv_formula, data = test)[ ,-1]  #drop the intercept
cox.train <- coxph(surv_formula, data = train)

cox_rp = as.matrix(x_test) %*% as.vector(unlist(coef(cox.train)))
cox_rp_train = as.matrix(x) %*% as.vector(unlist(coef(cox.train)))
#-------------------
metrics_list <- get_survmetrics(train = train, test = test, cox_rp = cox_rp, 
                                cox_rp_train = cox_rp_train)
models <- c(models, "Cox's-genomic-cluster")
cindex <- c(cindex, metrics_list[["Cindex.cox"]])
bs <- c(bs, metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(auc, metrics_list[["AUC.cox"]]$AUC)
D <- c(D, metrics_list[["D.cox"]]$d.index)
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc)
## ClinicoGenomic models
surv_formula <- as.formula(paste0("Surv(time, status) ~ ", paste(c(clinical_vars, "cluster") , collapse = "+") ) )

## Create model matrix --
x <- model.matrix(surv_formula, data = train)[ ,-1]  #drop the intercept
x_test <- model.matrix(surv_formula, data = test)[ ,-1]  #drop the intercept
cox.train <- coxph(surv_formula, data = train)

cox_rp = as.matrix(x_test) %*% as.vector(unlist(coef(cox.train)))
cox_rp_train = as.matrix(x) %*% as.vector(unlist(coef(cox.train)))
metrics_list <- get_survmetrics(train = train, test = test, cox_rp = cox_rp, 
                                cox_rp_train = cox_rp_train)
#-------------------
models <- c(models, "Cox's-clinico-cluster-genomic")
cindex <- c(cindex, metrics_list[["Cindex.cox"]])
bs <- c(bs, metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(auc, metrics_list[["AUC.cox"]]$AUC)
D <- c(D, metrics_list[["D.cox"]]$d.index)
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc, Dindex = D)

out <- list("metrics_test" = metrics_test, "k" = k)

saveRDS(out, paste0("~/iClusterPlusUtils/lusc/results/metrics_cluster", i, ".RDS") )
