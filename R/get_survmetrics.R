# Function to predict the survival for each individual
#'
#' Predict survival from a stan_surv object
#'
#' @export get_metrics
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
  bs.cox <- survcomp::sbrier.score2proba(data.tr = data.frame(time = train$time, event = train$status, score = cox_rp_train), data.ts = data.frame(time = test$time, event = test$status, score = cox_rp), method = c("prodlim"))

  
  #AUC ROC
  AUC.cox <- survcomp::tdrocc(x = as.vector(cox_rp) , surv.time = test$time, surv.event = test$status,  time = max(obs.test.time) )
  
  list("Cindex.cox" = Cindex.cox,  "bs.cox" = bs.cox, "AUC.cox" = AUC.cox)
}

