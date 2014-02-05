print.mcrot <-
function(x, digits= max(3, getOption("digits") - 3), ...){  
  cat("\t mcroast: Multiple Contrast Rotation Tests\n\n")
  v <- sqrt(diag(x$K %*% solve(t(x$X) %*% x$X) %*% t(x$K)))
  rnames <- paste(rep(rownames(x$B), each=ncol(x$B)), rep(colnames(x$B), times=nrow(x$B)), sep=" @ ")
  dd <- data.frame(Estimates=as.vector(t(x$B)),
                   Std.Err=as.vector(t(x$SD*v)),
                   Statistic=as.vector(t(x$statistic)),
                   "p-values"=as.vector(t(x$pvalue)),
                   row.names=rnames)
  printCoefmat(dd, digits=digits, cs.ind=1:2, tst.int=3, has.Pvalue=TRUE)
}

