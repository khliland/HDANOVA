contr.weighted <- function (x, base){
  if(missing(base))
    base <- levels(x)[1]
  frequencies <- table(x)
  base <- which(levels(x) == base)
  contr <- contr.treatment(length(frequencies), base = base)
  contr[base, ] <- -1 * frequencies[-base]/frequencies[base]
  dimnames(contr) <- list(names(frequencies),names(frequencies[-base]))
  return(contr)
}
#   dv <- dummyvar(x)
#   contr <- crossprod(dv*rep(median(sqrt(colSums(dv)))/colSums(dv), each=nrow(dv)))
#   contr[base, ] <- -contr[base,base]
#   browser()
#   contr[,-base]
# }
