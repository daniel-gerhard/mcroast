confint.mcrot <- function(object, parm, level = 0.95, ...){
  K <- object$K
  X <- object$X
  v <- sqrt(diag(K %*% solve(t(X) %*% X) %*% t(K)))
  quant <- quantile(object$mrot, level)
  upp <- object$B + quant*object$SD * v 
  low <- object$B - quant*object$SD * v
  
  out <- list()
  out$B <- object$B
  out$level <- level
  out$lower <- low
  out$upper <- upp
  class(out) <- "mcroastci"
  return(out)
}