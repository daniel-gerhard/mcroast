mixedroastlmerlist <- function(lmerlist, K, nrot=9999, adjusted=TRUE){
  Vl <- lapply(lmerlist, function(m) getME(m, "Z") %*% (sigma(m)^2 * getME(m, "Lambda") %*% getME(m, "Lambdat")) %*% getME(m, "Zt") + diag(nrow(model.matrix(m)))*sigma(m)^2) 
  Y <- sapply(1:length(lmerlist), function(i) as.vector(solve(t(chol(Vl[[i]]))) %*% getME(lmerlist[[i]], "y")))
  y <- t(Y)
  ng <- nrow(y)
  n <- ncol(y)
  index <- rep.int(TRUE, ng) 
  X <- model.matrix(lmerlist[[1]])
  design <- as.matrix(X)
  p <- ncol(design)
  p0 <- p - 1
  d <- n - p
  
  namey <- names(lmerlist)
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
    Xl <- lapply(Vl, function(V) contrastAsCoef(as.matrix(solve(t(chol(V))) %*% design), contrast, first = FALSE)$design)
    qrl <- lapply(Xl, function(Xc) qr(Xc))
    signc <- sapply(qrl, function(qr) sign(qr$qr[p, p]))
    effects <- sapply(1:nrow(y), function(i) qr.qty(qrl[[i]], t(y[i,,drop=FALSE])))
    s2 <- colMeans(effects[-(1:p), , drop = FALSE]^2)
    sd.post <- sqrt(s2)
    d0 <- 0
    
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
    sdr.post <- sqrt(s2r) 
    
    modtro <- signc * Br/sdr.post
    modtr <- zscoreT(modtro, df = d0 + d)
    return(list(modt=modt, modtr=modtr, B=signc * B, sdpost=sd.post))
  })
  
  Bmat <- t(sapply(outl, function(x) x$B) * t(sapply(Vl, function(V){
    Xv <- as.matrix(solve(t(chol(V))) %*% X)
    sqrt(diag(K %*% solve(t(Xv) %*% Xv) %*% t(K)))
  })))
  
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
  
  stderr <- sapply(lmerlist, function(mm) sqrt(diag(K %*% as.matrix(vcov(mm)) %*% t(K))))
  
  colnames(Bmat) <- colnames(pv) <- rownames(modt) <- colnames(stderr) <- namey
  rownames(Bmat) <- rownames(pv) <- colnames(modt) <- rownames(stderr) <- namek  
  
  
  
  out <- list()
  out$Y <- t(y)
  out$X <- model.matrix(lmerlist[[1]])
  out$statistic <- t(modt)
  out$pvalue <- pv 
  out$B <- Bmat
  out$mrot <- mrot
  out$K <- K
  out$stderr <- stderr
  
  class(out) <- "mixedrot"
  return(out)
}