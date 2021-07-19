#' Loglikelihood
#'
#' Compute the loglikelihood of a CLCM
#' @param X matrix of item responses
#' @param item.type character vector of item types
#' @param eta matrix operationalizing the condensation rules
#' @param categories.j numeric vector of number of categories per item
#' @return matrix of N by 2^K
#' @export

compute_loglike <- function(X, item.type, param, eta, categories.j){

  loglike <- matrix(0, nrow = nrow(X), ncol = nrow(eta))

  for(j in 1:length(item.type)){

    loglike_j <- matrix(0, nrow = nrow(X), ncol = 2)
    param_j <- param[[j]]
    item.type_j <- item.type[j]
    categories.j_j <- categories.j[[j]]

    logPj <- Return_logPj(item.type_j = item.type_j, param_j = param_j, j = j, eta = eta, categories.j_j = categories.j_j)

    if(item.type_j == 'Normal'){

      mean.l <- logPj[grep('mean', names(logPj))] # means of the latent classes
      sigma <- logPj['sigma'] # sd pooled across latent classes
        for(i in 1:nrow(loglike_j)){
            ##loglike_j[i, ] <- pnorm(X[i, j], mean = mean.l, sd = sigma, lower.tail = TRUE)
            loglike_j[i, ] <- dnorm(X[i, j], mean = mean.l, sd = sigma)
        }# end i loop

      loglike_j <- log(loglike_j)


    } else {

      if(item.type_j == 'Beta'){
        shape1.l <- logPj[grep('shape1', names(logPj))] # means of the latent classes
        shape2 <- logPj['shape2'] # shape2 - not differentiated across latent classes - legit??
          for(i in 1:nrow(loglike_j)){
              ## loglike_j[i, ] <- pbeta(X[i, j], shape1 = shape1.l, shape2 = shape2, lower.tail = TRUE)
              loglike_j[i, ] <- dbeta(X[i, j], shape1 = shape1.l, shape2 = shape2)
          }# end i loop

        loglike_j <- log(loglike_j)

      } else {

        # Categorical Items - Ordinal, Count, etc.
        for(i in 1:nrow(loglike_j)){
            loglike_j[i, ] <- logPj[ X[i, j] + 1,  ]
            }# end i loop

    }}# end if else - expand with different item types


    # loglike_j will be Nx2, because conjunctive condensation rule reduces 2^K down to 2
    # Expand back out to 2^K
    loglike_j <- loglike_j %*% rbind(1-eta[, j], eta[, j])
    # Now loglike_j is Nx2^K

    loglike_j <- replace(loglike_j, is.na(loglike_j), 0)
    loglike <- loglike + loglike_j


  }#end j loop

  return(loglike)

}#end compute_loglike
