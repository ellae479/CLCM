
compute_post <- function(X, item.type, param, lprior, eta, categories.j){

  loglike <- compute_loglike(X, item.type, param, eta, categories.j)

  # Normalize posterior distributions
  post.updated <- loglike + lprior
  post.updated <- exp(post.updated)
  post.updated <- post.updated/matrix(rowSums(post.updated, na.rm = T),
                                      nrow = nrow(post.updated),
                                      ncol = ncol(post.updated), byrow = F) #normalize

  return(post.updated)

}#end compute_post
