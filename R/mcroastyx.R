mcroastyx <-
function(Y, X, K, nrot=9999, adjusted=TRUE, moderated=FALSE){
  y <- as.matrix(t(Y))
  ng <- nrow(y)
  n <- ncol(y)
  index <- rep.int(TRUE, ng)
  design <- as.matrix(X)
  p <- ncol(design)
  p0 <- p - 1
  d <- n - p
  
  namey <- colnames(Y)
  namek <- rownames(K)
  if (is.null(namey)) namey <- paste("Y", 1:ncol(Y), sep="")
  if (is.null(namek)) namek <- paste("C", 1:nrow(K), sep="") 
  
  if (adjusted){
    crm <- cov2cor(K %*% qr.solve(crossprod(X)) %*% t(K))
    ssvd <- svd(crm)
    rval <- t(ssvd$v %*% (t(ssvd$u) * sqrt(ssvd$d)))
    R <- t(replicate(nrot, as.vector(matrix(rnorm((d+1) * ncol(crm)), nrow = (d+1), byrow = TRUE) %*% rval)))    
    Rid <- matrix(as.logical(rep(1, d+1) %o% diag(nrow(K))), ncol=nrow(K))
  }
  
  outl <- lapply(1:nrow(K), function(ir){
    contrast <- K[ir,]    
    X <- contrastAsCoef(design, contrast, first = FALSE)$design
    qr <- qr(X)
    signc <- sign(qr$qr[p, p])
    effects <- qr.qty(qr, t(y))
    s2 <- colMeans(effects[-(1:p), , drop = FALSE]^2)
    if (moderated == TRUE){
      sv <- squeezeVar(s2, df = d)
      d0 <- sv$df.prior
      s02 <- sv$var.prior
      effects <- effects[, index, drop = FALSE]
      sd.post <- sqrt(sv$var.post[index])
    } else {
      sd.post <- sqrt(s2)
      d0 <- 0
    }
    
    nset <- ncol(effects)
    if (p0 > 0) Y <- effects[-(1:p0), , drop = FALSE] else Y <- effects
    YY <- colSums(Y^2)
    B <- Y[1, ]
    modto <- signc * B/sd.post
    modt <- zscoreT(modto, df = d0 + d)
    
    if (adjusted){
      Rr <- R[, Rid[,ir]]      
    } else {
      Rr <- matrix(rnorm(nrot * (d + 1)), nrot, d + 1)
    }
    Rr <- Rr/sqrt(rowSums(Rr^2))
    Br <- Rr %*% Y
    s2r <- (matrix(YY, nrot, nset, byrow = TRUE) - Br^2)/d
    if (moderated == TRUE){
      if (is.finite(d0)){
        sdr.post <- sqrt((d0 * s02 + d * s2r)/(d0 + d)) 
      } else {
        sdr.post <- sqrt(s02)
      }
    } else {
      sdr.post <- sqrt(s2r) 
    }
    modtro <- signc * Br/sdr.post
    modtr <- zscoreT(modtro, df = d0 + d)
    return(list(modt=modt, modtr=modtr, B=signc * B, sdpost=sd.post))
  })
  
  Bmat <- t(sapply(outl, function(x) x$B)) * sqrt(diag(K %*% solve(t(X) %*% X) %*% t(K)))
  SDmat <- t(sapply(outl, function(x) x$sdpost))
  
  if (adjusted){
    modt <- sapply(outl, function(x) x$modt)
    modtr <- matrix(simplify2array(lapply(outl, function(x) x$modtr)), nrow=nrot)
    mrot <- apply(abs(modtr), 1, max)
    tcomp <- sapply(mrot, function(x) x >= as.vector(t(abs(modt))))
    pv <- matrix((rowSums(tcomp) + 1)/(nrot + 1), nrow=nrow(K))
  } else {
    mrot <- sapply(outl, function(x) apply(abs(x$modtr), 1, max))
    modt <- sapply(outl, function(x) abs(x$modt))    
    rsum <- sapply(1:ncol(modt), function(i) rowSums(sapply(mrot[,i], function(x) x >= modt[,i])))
    pv <- t((rsum + 1)/(nrot + 1))
  }    
  colnames(Bmat) <- colnames(SDmat) <- colnames(pv) <- rownames(modt) <- namey
  rownames(Bmat) <- rownames(SDmat) <- rownames(pv) <- colnames(modt) <- namek  
  
  out <- list()
  out$Y <- t(y)
  out$X <- X
  out$statistic <- t(modt)
  out$pvalue <- pv 
  out$B <- Bmat
  out$SD <- SDmat
  out$mrot <- mrot
  out$K <- K
  
  class(out) <- "mcrot"
  return(out)
}
