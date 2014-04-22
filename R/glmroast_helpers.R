extractglmweights <- function(object){
  family <- family(object)
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta  
  mf <- model.frame(object)
  est <- coefficients(object)
  X <- model.matrix(object)
  nobs <- nrow(X)
  weights <- model.weights(object)
  if (is.null(weights)) weights <- rep.int(1, nobs)
  y <- model.response(mf)    
  eta <- X %*% est
  etastart <- eta
  eval(family$initialize)  
  mu <- linkinv(eta)
  mu.eta.val <- mu.eta(eta)  
  offset <- model.offset(mf)
  if (!length(offset)) offset <- rep(0, nobs)
  as.vector(sqrt((weights * mu.eta.val^2)/variance(mu)))
}

glmresponse <- function(object){
  family <- family(object)
  variance <- family$variance
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta  
  mf <- model.frame(object)
  est <- coefficients(object)
  X <- model.matrix(object)
  nobs <- nrow(X)
  weights <- model.weights(object)
  if (is.null(weights)) weights <- rep.int(1, nobs)
  y <- model.response(mf)  
  eta <- X %*% est
  etastart <- eta
  eval(family$initialize)
  mu <- linkinv(eta)
  mu.eta.val <- mu.eta(eta)  
  offset <- model.offset(mf)
  if (!length(offset)) offset <- rep(0, nobs)
  z <- (eta - offset) + (y - mu)/mu.eta.val
  w <- as.vector(sqrt((weights * mu.eta.val^2)/variance(mu)))
  return(z*w)
}
