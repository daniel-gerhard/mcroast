mixedroast <- function(formula, data=NULL, K, nrot = 9999, adjusted = TRUE){
  ymat <- model.response(model.frame(as.formula(paste(deparse(formula[[2]], width.cutoff=500), "~ 1")), data=data))
  ynam <- colnames(ymat)
  forms <- sapply(ynam, function(yna) as.formula(paste(yna, "~", paste(deparse(formula[[3]]), collapse=""))))
  lmerlist <- lapply(forms, function(ff) lmer(ff, data=data))
  mixedroastlmerlist(lmerlist=lmerlist, K=K, nrot=nrot, adjusted=adjusted)
}