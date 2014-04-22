print.glmrot <- function(x, digits= max(3, getOption("digits") - 3), ...){  
    cat("\t mcroast: Multiple Contrast Rotation Tests\n\n")
    rnames <- paste(rep(rownames(x$B), each=ncol(x$B)), rep(colnames(x$B), times=nrow(x$B)), sep=" @ ")
    dd <- data.frame(Estimates=as.vector(t(x$B)),
                     Statistic=as.vector(t(x$statistic)),
                     "p-values"=as.vector(t(x$pvalue)),
                     row.names=rnames)
    printCoefmat(dd, digits=digits, cs.ind=1, tst.int=2, has.Pvalue=TRUE)
  }