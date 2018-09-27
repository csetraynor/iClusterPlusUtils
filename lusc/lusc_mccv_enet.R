#mc-cross-validation.R
# Load libraries
suppressMessages(library(dplyr))
# For glmnet
suppressMessages(library(glmnet))
suppressMessages(library(glmnetUtils))

# For survival analysis
library(survival)
library(survcomp)
# For MC-crossvalidations
suppressMessages(library(rsample))
suppressMessages(library(tidyposterior))

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

## Add patient id
gene$matrix <- dplyr::as_data_frame(gene$matrix, rownames = "sample_id")
meth450$matrix <- dplyr::as_data_frame(meth450$matrix, rownames = "sample_id")
seg.cn <- dplyr::as_data_frame(seg.cn, rownames = "sample_id")
## join genomic head vars
genomic <- dplyr::left_join(gene$matrix, meth450$matrix, by = "sample_id") %>% 
  dplyr::left_join(., seg.cn, by = "sample_id")
rm(gene)
rm(meth450)
rm(seg.cn)
genomic_true_vars <- colnames(genomic %>%
                           dplyr::select(-"sample_id"))
genomic_vars <- paste0("var", seq_along(1:length(genomic_true_vars)))
colnames(genomic)[-match("sample_id", colnames(genomic))] <- genomic_vars
# get clinical data -----
lusc_mccv <- readRDS("~/iClusterPlusUtils/data/mc_lusc_clinical.RDS")
vfold <- lusc_mccv$splits[[i]] #choose splits
train_id <- vfold$in_id

train <- rsample::analysis(vfold)
test <- rsample::assessment(vfold)


## Join genomic and clinical data
train <- dplyr::left_join(train, sample %>% dplyr::select(patient_id, sample_id), by = "patient_id") %>% dplyr::left_join(., genomic,  by = "sample_id") %>% dplyr::select(-sample_id)

test <- dplyr::left_join(test, sample %>% dplyr::select(patient_id, sample_id), by = "patient_id") %>% dplyr::left_join(., genomic,  by = "sample_id") %>% dplyr::select(-sample_id)

nzvars <- caret::nearZeroVar(train[ ,genomic_vars])
if (length(nzvars) > 0) {
  train[ ,genomic_vars] <- train[ ,genomic_vars][ ,-nzvars]
  genomic_vars <- genomic_vars[-nzvars]
}
test <- test[ ,colnames(train)]

### Standardise covariates ----------------------
preProcValues <- caret::preProcess(train[ ,c("age", genomic_vars) ], method = c("center", "scale") )
trainTransformed <- as.data.frame( predict(preProcValues, train[c("age", genomic_vars) ] ) )
colnames(trainTransformed) <- c("age", genomic_vars)
#colnames(trainTransformed) = c("age_std")
train <- cbind(train[ ,-match(c("age", genomic_vars), colnames(train) )], trainTransformed) %>%
  dplyr::mutate(time = os_months,
                status = os_status)

rm(trainTransformed) # but we keep preProcValues
## transform test set with training vars
testTransformed <- as.data.frame( predict(preProcValues, test[ ,c("age", genomic_vars)]) )
test <- cbind(test[ ,-match(c("age", genomic_vars), colnames(test) )], testTransformed) %>%
  dplyr::mutate(time = os_months,
                status = os_status)
# delete 
rm(list = c("testTransformed"))

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
  
  #AUC ROC
  AUC.cox <- survcomp::tdrocc(x = as.vector(cox_rp) , surv.time = test$time, surv.event = test$status,  time = max(obs.test.time) )
  
  list("Cindex.cox" = Cindex.cox,  "bs.cox" = bs.cox, "AUC.cox" = AUC.cox)
}

metrics_list <- get_survmetrics(train = train, test = test, cox_rp = cox_rp, 
                            cox_rp_train = cox_rp_train)
cindex <- c(metrics_list[["Cindex.cox"]])
bs <- c(metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(metrics_list[["AUC.cox"]]$AUC)
models <- c("Cox's-clinical")
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc)
## Genomic models
x <- train[, genomic_vars]
x_test <- test[, genomic_vars]
y <- train[, c("time", "status")]
y$status <- as.integer(y$status)
y <- as.matrix(y)
x <- as.matrix(x)
y_test <- test[, c("time", "status")]
in_data <- cbind.data.frame(y, x)
in_data_test <- cbind.data.frame(y_test, x_test)
# expected number of non-zeros
nz = ceiling(0.01 * length(genomic_vars))
z = length(genomic_vars) - nz
# register for parallelisation if available
set.seed(82405)
foldid <- caret::createFolds(in_data$status, k = 10, list = FALSE)
## do crossvalidation to find optimal alpha and lambda
elasticnet_stuff <- cva.glmnet(x, y, family = "cox", 
                               grouped = TRUE, foldid = foldid, parallel = FALSE)
# calculate prognostic index
cvm <- sapply(elasticnet_stuff$modlist, function(i) min(i$cvm))
elasticnet <- elasticnet_stuff$modlist[[match(min(cvm), cvm)]]
lam <- match(min(cvm), elasticnet$cvm)

# calculate risk score
enet_rp <- as.matrix(x_test) %*% as.vector(coef(elasticnet$glmnet.fit)[, 
                                                                       lam])
enet_rp_train <- as.matrix(x) %*% as.vector(coef(elasticnet$glmnet.fit)[, 
                                                                        lam])
#-------------------
metrics_list <- get_survmetrics(train = train, test = test, cox_rp = enet_rp, 
                                cox_rp_train = enet_rp_train)
models <- c(models, "Cox's-enet-genomic")
cindex <- c(cindex, metrics_list[["Cindex.cox"]])
bs <- c(bs, metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(auc, metrics_list[["AUC.cox"]]$AUC)
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc)

# ## Clinico-genomic model give dummy names to vars to avoid
# issues when creating model formula
x_genomic <- train[, genomic_vars]
x_genomic_test <- test[, genomic_vars]
## prepare input data
x <- cbind(x_clinical, x_genomic)
x_test <- cbind(x_clinical_test, x_genomic_test)
y <- train[, c("time", "status")]
y_test <- test[, c("time", "status")]
in_data <- cbind.data.frame(y, x)
## Train model register for parallelisation if available
set.seed(1347294)
foldid <- caret::createFolds(train$status, k = 10, list = FALSE)
pf = rep(1, ncol(x))
pf[match(colnames(x_clinical), colnames(x))] <- 0  #clinico-genomic assumption
# do crossvalidation to find optimal alpha and lambda
elasticnet_stuff <- cva.glmnet(as.matrix(x), as.matrix(y), family = "cox", 
                               grouped = TRUE, foldid = foldid, parallel = FALSE, penalty.factor = pf)
## 
cvm <- sapply(elasticnet_stuff$modlist, function(i) min(i$cvm))
elasticnet <- elasticnet_stuff$modlist[[match(min(cvm), cvm)]]
lam <- match(min(cvm), elasticnet$cvm)
## calculate risk score
enet_rp_train <- as.matrix(x) %*% as.vector(coef(elasticnet$glmnet.fit)[, 
                                                                        lam])
enet_rp <- as.matrix(x_test) %*% as.vector(coef(elasticnet$glmnet.fit)[, 
                                                                       lam])
in_data_test <- cbind(y_test, x_test)
metrics_list <- get_survmetrics(train = in_data, test = in_data_test, 
                            cox_rp = enet_rp,
                            cox_rp_train = enet_rp_train)
#-------------------
models <- c(models, "Cox's elastic net: clinico-genomic")
cindex <- c(cindex, metrics_list[["Cindex.cox"]])
bs <- c(bs, metrics_list[["bs.cox"]]$bsc.integrated)
auc <- c(auc, metrics_list[["AUC.cox"]]$AUC)
metrics_test <- data.frame(models = models, cindex = cindex, bs = bs, 
                           auc = auc)


saveRDS(metrics, paste0("~/iClusterPlusUtils/lusc/results/metrics_cluster", i, ".RDS") )
