adjustHessian <- function(Hessian){
#### make sure matrix Hessian is inversible
    rank.m = rankMatrix(Hessian)
    if(rank.m < ncol(Hessian)){
      tmp.svd = svd(Hessian)
      Hessian = (tmp.svd$u%*%diag(tmp.svd$d)%*%t(tmp.svd$v) -
                   0.01*tmp.svd$u%*%diag(c(rep(0, rank.m),
                                           rep(1, ncol(Hessian) - rank.m))
                                         )%*%t(tmp.svd$v)
                 )
      }
    return(Hessian)
}
