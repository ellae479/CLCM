

compute_lprior <- function(post, K, Z = NULL, type = NULL, reg.formula = NULL){


  post <- as.matrix(post)
  # Do all this shit later:
  # if(!is.null(Z)){ #if you passed covariates:
  #   if(nrow(Z) != nrow(post)){stop('Check your sample size of covariates and posterior distribution')}
  #     if(is.null(type)){stop('Specify the prior distribution you want to estimate - see doc')}else{ #if you passed covariates AND you passed a 'type', then...
  #       if(type == 'max_LL' & is.null(loglike)){ stop('Need to pass the loglike for this estimation approach')}
  #     }}
  #


########################################################
  # Different types of prior estimation:

  # No Covariate, standard empirical Bayes procedure:
  if(type == 'standard_EB'){

      lprior <- matrix(log(colMeans(post)), nrow = nrow(post), ncol = 2^K, byrow = T)

  } # end "standard_EB"



  # Empirical Bayes - implement groups by taking mean of the posterior distribution
  ## similar to standard EB approach but stratified on covariate subgroups

  if(type == 'groups_EB'){

      lprior <- matrix(NA, nrow = nrow(post), ncol = ncol(post))
      X <- model.matrix(reg.formula, data = as.data.frame(Z))
      groups <- apply(X, 1, paste0, collapse = '')

          for(i in sort(unique(groups)) ){

            tmp <- colMeans(post[ which(groups == i), ], na.rm = T)
            lprior[ which(groups == i), ] <- matrix(tmp, nrow = sum(groups == i), ncol = ncol(lprior), byrow=T)

          }

  lprior <- log(lprior)

  } # end "groups_EB"



 # Note that max_LL can be used with any type of covariate
  if(type == 'latent_regression'){


    XX <- model.matrix(as.formula(paste0(' ~ ', reg.formula)), data = as.data.frame(Z))
    vP <- rep(0, ncol(XX)*(2^K-1))

    # to estimate using this approach, first load the function:
       full_LL <- function(vP, post, XX){

              #print(vP)
              vP <- matrix(vP, nrow = ncol(XX), ncol = 2^K-1, byrow = F)
              XB <- XX %*% vP
              p <- exp(XB)/(1 + apply(exp(XB), 1, sum))
              p <-  cbind(1 - rowSums(p), p)
              LL <- post * log(p)
              ##LL <- sum(log(rowSums(exp(loglike) * p, na.rm = T)), na.rm = T) - notice how this is different
              # original was optimizing the actual objective function
              # not sure if it's correct, so leave it here
              LL <- -1*LL # optimization will minimize function
              LL <- sum(LL)
              return(LL)
          }#end function full_LL()

    out <- optim(par = vP, fn = full_LL, gr = NULL, method = "BFGS", post = post, XX = XX)

      vP <- matrix(out$par, nrow = ncol(XX), ncol = 2^K-1, byrow = F, dimnames = list(colnames(XX), NULL))
      XB <- XX %*% vP
      p <- exp(XB)/(1 + apply(exp(XB), 1, sum))
      p <-  cbind(1 - rowSums(p), p)
      lprior <- log(p)
      lprior <- list('vP' = vP, 'lprior' = lprior) #


  }#end type == 'latent_regression'


   # Note that max_LL can be used with any type of covariate
  if(type == 'multinom_latent_regression'){


    XX <- model.matrix(as.formula(paste0(' ~ ', reg.formula)), data = as.data.frame(Z))
    vP <- rep(0, ncol(XX)*(2^K-1)) # actually I think this was always correct

    # to estimate using this approach, first load the function:
       full_LL <- function(vP, post, XX){

              #print(vP)
              vP <- matrix(vP, nrow = ncol(XX), ncol = 2^K-1, byrow = F)
              XB <- XX %*% vP
              p <- exp(XB)/(1 + apply(exp(XB), 1, sum))
              p <-  cbind(1 - rowSums(p), p)
              LL <- post * log(p)
              LL <- -1*LL # optimization will minimize function
              LL <- sum(LL)
              return(LL)
          }#end function full_LL()

    out <- optim(par = vP, fn = full_LL, gr = NULL, method = "BFGS", post = post, XX = XX)

      vP <- matrix(out$par, nrow = ncol(XX), ncol = 2^K-1, byrow = F, dimnames = list(colnames(XX), NULL))
      XB <- XX %*% vP
      p <- exp(XB)/(1 + apply(exp(XB), 1, sum))
      p <-  cbind(1 - rowSums(p), p)
      lprior <- log(p)
      lprior <- list('vP' = vP, 'lprior' = lprior) #


  }#end type == 'multinom_latent_regression'


  return(lprior)

}# end lprior function
