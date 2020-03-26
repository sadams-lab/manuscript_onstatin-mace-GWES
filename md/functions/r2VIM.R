# Based on https://github.com/silkeszy/Pomona

var.sel.r2vim <- function(formula = formula, data = data, no.runs = 10, factor = 1, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                          no.threads = 1, method = "ranger", type = "regression", case.weights = NULL) {
  
  ## importance for each run 
  imp.all = NULL
  for (r in 1:no.runs) {
    print(paste("run", r))
    rf = wrapper.rf(formula = formula, data = data,
                    ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop, no.threads = no.threads,
                    method = method, type = type, case.weights = case.weights)
    imp.all = cbind(imp.all, get.vim(rf))
  }
  
  ## factors
  min.global = min(imp.all)
  if (min.global >= 0) {
    stop("Global minimal importance score is not negative!")
  }
  no.neg.min = 0
  fac = matrix(nrow = nrow(imp.all), ncol = ncol(imp.all),
               dimnames = dimnames(imp.all))
  for (i in 1:ncol(imp.all)) {
    x = imp.all[,i]
    min = min(x)
    if (min >= 0) {
      no.neg.min = no.neg.min + 1
      fac[,i] = x / abs(min.global)
    } else {
      fac[, i] = x / abs(min)
    }
  }
  if (no.neg.min > 0) {
    print(paste(no.neg.min, "runs with no negative importance score!"))
  }
  fac.min = apply(fac, 1, min)
  fac.med = apply(fac, 1, median)
  
  ## select variables
  ind.sel = as.numeric(fac.min >= factor)
  
  ## info about variables
  info = data.frame(imp.all, fac, fac.min, fac.med, ind.sel)
  colnames(info) = c(paste("vim.run.", 1:no.runs, sep = ""),
                     paste("rel.vim.run.", 1:no.runs, sep = ""),
                     "rel.vim.min", "rel.vim.median", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

wrapper.rf <- function(formula = formula, data = data, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       method = "ranger", type = "regression", case.weights = case.weights, ...) {
  
  
  ## set global parameters
  nodesize = floor(nodesize.prop * nrow(data@gtdata))
  mtry = floor(mtry.prop * ncol(data@gtdata))
  if (mtry == 0) mtry = 1
  
  
  
  rf = ranger::ranger(formula = formula, 
                      data = data,
                      probability = FALSE,
                      importance = "permutation", 
                      scale.permutation.importance = FALSE,
                      num.trees = ntree,
                      mtry = mtry,
                      min.node.size = nodesize,
                      num.threads = no.threads,
                      write.forest = TRUE,
                      case.weights = case.weights,
                      ...)
  return(rf)
}


calculate.error <- function(rf, true, test.set = NULL) {
  
  if (is(rf, "ranger")) {
    if (!is.null(test.set)) {
      pred = predict(rf, data = test.set)$predictions
    } else {
      pred = rf$predictions
    }
    if (rf$treetype == "Probability estimation") {
      pred = pred[, 2]
    }
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  
  if ((is(rf, "randomForest") && rf$type == "classification") |
      (is(rf, "ranger") && rf$treetype == "Classification")) {
    conf.matrix = table(pred = pred, true = true)
    tp = conf.matrix[2, 2]
    tn = conf.matrix[1, 1]
    fn = conf.matrix[2, 1]
    fp = conf.matrix[1, 2]
    
    ## accuracy
    acc = (tp + tn) / sum(conf.matrix)
    
    ## Matthews correlation coefficient
    mcc = (tp * tn - fp * fn) /
      sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn))
    
    ## sensitivity
    sens = tp / (tp + fn)
    
    ## specificity
    spec = tn / (fp + tn)
    
    error = c(err = 1 - acc, acc = acc, mcc = mcc, sens = sens, spec = spec)
  } else {
    mse = sum((pred - true)^2, na.rm = TRUE) / sum(!is.na(pred))
    
    ## pseudo R-squared uses sum of squared differences divided by n instead of variance!
    v = sum((true - mean(true))^2) / length(true)
    rsq = 1 - mse/v
    error = c(rmse = sqrt(mse), rsq = rsq)
  }
  
  return(error)
}

get.vim <- function(rf) {
  if (is(rf, "ranger")) {
    vim = ranger::importance(rf)
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  return(vim)
}

