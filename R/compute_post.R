#' Posterior distribution
#'
#' Compute the posterior distribution of a CLCM
#' @param X matrix of item responses
#' @param item.type character vector of item types
#' @param param list of item parameters
#' @param lprior N by 2^K matrix of log of prior distribution
#' @param eta matrix operationalizing the condensation rules
#' @param categories.j numeric vector of number of categories per item
#' @return matrix of N by 2^K
#' @export
#'

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
