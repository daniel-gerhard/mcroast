confint.mixedrot <- function(object, parm, level = 0.95, ...){
  K <- object$K
  X <- object$X
  v <- object$stderr
  quant <- quantile(object$mrot, level)
  upp <- object$B + quant * v 
  low <- object$B - quant * v
  
  out <- list()
  out$B <- object$B
  out$level <- level
  out$lower <- low
  out$upper <- upp
  class(out) <- "mcroastci"
  return(out)
}