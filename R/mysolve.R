mysolve <- function(A){
  ### find the inverse of matrix A
  if(det(A) < 1e-10){
    tmp.svd = svd(A)
    A = (tmp.svd$u%*%diag(tmp.svd$d)%*%t(tmp.svd$v) -
           0.1*tmp.svd$u%*%diag(
             c(rep(0, sum(tmp.svd$d > 1e-15)),
               rep(1, sum(tmp.svd$d < 1e-10))))%*%t(tmp.svd$v)
    )
  }
  solve(A)
}