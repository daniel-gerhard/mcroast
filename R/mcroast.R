mcroast <-
function(formula, data=NULL, K, nrot=9999, adjusted=TRUE, moderated=FALSE){    
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")  
  Y <- model.response(mf, "numeric")
  X <- model.matrix(mt, mf)    
  mc <- mcroastyx(Y, X, K=K, nrot=nrot, adjusted=adjusted)
}
