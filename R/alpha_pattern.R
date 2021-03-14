
pattern<-function(K){
  alpha=vector()
  for (k in K:1){
    tmp=vector()
    base=rbind( matrix(0, nrow=2^(K-k),ncol=1), matrix(1,nrow=2^(K-k),ncol=1) )
    base=rep(base,2^(k-1))
    alpha=cbind(alpha,base)}

  alpha=cbind(alpha, rowSums(alpha))
  alpha=alpha[(order(alpha[,(K+1)])),]
  alpha=alpha[,-(K+1)]
  alpha<-as.matrix(alpha)
  return(alpha)
}#end patttern


# DINA 'eta' - creates DINA model:
# alpha<-pattern(K)
# eta=alpha %*% t(Q)
# eta=ifelse(eta == matrix(1,2^K,1) %*% colSums(t(Q)),1, 0 )
