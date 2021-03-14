
compute_like_jacobian <- function(param, X, item.type, eta, categories.j, post, item.names){

  # How to take unlisted parameters and re-create the full list??
    vP <- param
    param.vP <- vector(mode="list", length =length(item.type))

  for(j in 1:length(item.type)){

      vP.j <- vP[ grep(item.names[j], names(vP), value = T) ]
      names(vP.j) <- sub(pattern = paste0(item.names[j], '.'), replacement = '', x = names(vP.j))
      param.vP[[j]] <- vP.j

  }

    names(param.vP) <- item.names
    #all(unlist(param.vP) == unlist(item.param))  # CHECK


  # Compute loglikelihood
  loglike <- compute_loglike(X = X,
                             item.type = item.type,
                             param = param.vP,
                             eta = eta,
                             categories.j = categories.j)

  like <- exp(loglike)
  likeNby1 <- like %*% colMeans(post) # P(X) - probability of observing that item response pattern given the parameter estimates
  return(likeNby1)

} #end function



