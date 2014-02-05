print.mcroastci <- function(x, digits= max(3, getOption("digits") - 3), ...){
  cat(paste("\t ", x$level*100, "% Confidence Intervals based on Multiple Contrast Rotation Tests\n\n", sep=""))
  rnames <- paste(rep(rownames(x$B), each=ncol(x$B)), rep(colnames(x$B), times=nrow(x$B)), sep=" @ ")
  dd <- data.frame(Estimates=as.vector(t(x$B)),
                   Lower=as.vector(t(x$lower)),
                   Upper=as.vector(t(x$upper)),
                   row.names=rnames)
  print(dd, digits=digits)
}