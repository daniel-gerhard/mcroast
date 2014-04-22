glmroast  <- function(formula, data=NULL, K, family=gaussian(), nrot = 9999, adjusted = TRUE, responsenumber=NULL, ...){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")  
  Y <- model.response(mf)
  X <- model.matrix(mt, mf)   
  if (is.null(responsenumber)) responsenumber <- 1:ncol(Y)
  ny <- length(unique(responsenumber))
  if (length(responsenumber) > ncol(Y)) stop("More responsenumbers available than response vectors.")
  if (class(family)[1] == "family") family <- lapply(1:ny, function(i) family)
  if (class(family)[1] != "list" | length(family) != ny) stop("family needs to be of class family or a list of family object with length equal to the number of responses.")  
  
  gid <- unclass(factor(responsenumber))
  glmlist <- lapply(1:ny, function(i) glm(Y[,gid == i] ~ X-1, data=data, family=family[[i]], ...))  
  glmroastlist(glmlist=glmlist, K=K, nrot=nrot, adjusted=adjusted)
}
