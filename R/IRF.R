
Return_logPj <- function(item.type_j, param_j, j, eta, categories.j_j){

  # 2.17.21 - Updated for EM algorithm - returns item Probabilities, latent class by category
  # Use: passes Pj depending on which item.type you're using
  # Purpose: only write the IRF one time, here
  logPj <- NULL

    if(item.type_j == 'Ordinal'){       logPj <- Ordinal_IRF(param_j = param_j, j=j, eta = eta) }
    if(item.type_j == 'Neg_Binom'){     logPj <- Neg_Binom_IRF(param_j = param_j, j=j, eta = eta, categories.j_j = categories.j_j) }
    if(item.type_j == 'Poisson'){       logPj <- Poisson_IRF(param_j = param_j, j=j, eta = eta, categories.j_j = categories.j_j) }
    if(item.type_j == 'ZINB'){          logPj <- ZINB_IRF(param_j = param_j, j=j, eta = eta, categories.j_j = categories.j_j) }
    if(item.type_j == 'ZIP'){           logPj <- ZIP_IRF(param_j = param_j, j=j, eta = eta, categories.j_j = categories.j_j) }
    if(item.type_j == 'Normal'){        logPj <- Normal_IRF(param_j = param_j, j=j, eta = eta) }
    if(item.type_j == 'Beta'){          logPj <- Beta_IRF(param_j = param_j, j=j, eta = eta) }
    if(item.type_j == 'Nominal'){       logPj <- Multinomial_IRF(param_j = param_j, j=j, eta = eta) }
  # Notes:
  # Ordinal doesn't need categories.j because it can resolve that from the intercepts
  # Normal doesn't need categories.j because it doesn't have categories


  return(logPj) # C by 2^K

} # end function "Return_Pj"


#

# 3.8.21: Nominal items, why not?
Multinomial_IRF <- function(param_j, j, eta){

  # Optimization function passes a vector "vP", not a matrix - re-create the matrix here
  if(is.vector(param_j)){
    param_j <- matrix(param_j, nrow = 2) # row 1: intercepts, row 2: slopes
  }

# The matrix consolidates the parameters into the two reduced latent classes:
  XB <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = T) %*% param_j
  p <- exp(XB)/(1 + apply(exp(XB), 1, sum))
  p <- cbind(1 - rowSums(p), p) # rowSums(p) = 1
  rownames(p) <- paste0('reduced_LC_', 1:nrow(p))
  colnames(p) <- paste0('category_', 1:ncol(p))
  p <- t(p)
  logPj <- log(p)
  return(logPj)  # C by 2^K

}# end Multinomial_IRF()



Beta_IRF <- function(param_j,  j, eta, X){

  # This is so fucking stupid - why not just write the goddamn parameters this way in the first place??
  beta.param <- c(param_j['intercept'],
                  param_j['intercept'] + param_j['slope'],
                  param_j['shape2'])

  names(beta.param) <- c(paste0('shape1_', c(0,1) ), 'shape2') # 3.7.21
  return(beta.param)

}# end Beta_IRF()



Normal_IRF <- function(param_j,  j, eta, X){

  # Again, so dumb, just write it this way in the initialize item parameters, please
  norm.param <- c(param_j['intercept'],
                  param_j['intercept'] + param_j['slope'],
                  param_j['sigma'])
  names(norm.param) <- c(paste0('mean_', c(0,1) ), 'sigma') # 3.7.21: modified
  ## print(norm.param)
  return(norm.param)

}# end Normal_IRF()







# 2.17.21 - Updated for EM algorithm - returns item Probabilities, latent class by category
Ordinal_IRF <- function(param_j, j, eta){

  inter <- matrix(param_j[grep('inter', names(param_j))], nrow = 2, ncol = length(grep('inter', names(param_j))), byrow=T)
  slope <- matrix(c(0,1) * param_j['slope'], nrow = 2, ncol = ncol(inter), byrow = F)

  p <- slope + inter
  p <- exp(p)/(1+exp(p))
  p <- cbind(1, p, 0)

    Pj <- vector()
  for(cc in 1:(ncol(p) - 1) ) {
    Pj <- cbind(Pj,  p[ , cc] - p[ , cc+1] )
    }

    Pj <- t(Pj) # 2.20.21, Output dimensions should be: C by 2^K
    logPj <- log(Pj)

    return(logPj)  # C by 2^K

}# end Ordinal_IRF()



Neg_Binom_IRF <- function(param_j,  j, eta, X, categories.j_j){

  inter <- matrix(param_j[grep('inter', names(param_j))], nrow = 2, ncol = length(grep('inter', names(param_j))), byrow=T)
  slope <- matrix(c(0,1) * param_j['slope'], nrow = 2, ncol = ncol(inter), byrow = F)

  mu <- slope + inter
  mu <- exp(mu)
  size <- exp(param_j['size'])
  # 7.22.19, This is an easy way to constrain the size parameter
  # to be real & positive (see Bolker textbook pg 167 Ecological Models and Data in R)
  # Note that in the future you should then exponentiate the param_j['size'] because
  # that will be what the model actually used

  prob <- size/(size + mu)

  prob <- matrix(prob, nrow = length(categories.j_j), ncol = nrow(prob), byrow = T) # Number of categories by number of latent classes
  Xj <- matrix(categories.j_j, nrow = length(categories.j_j), ncol = ncol(prob), byrow = F) # same, C by 2^K
  # Note that the Xj is NOT the item responses, it's just the possible categories

  logPj <- lgamma(Xj + size) - lgamma(size) - lgamma(Xj + 1) + size*log(prob) + (Xj)*log(1-prob)

  # The probabilities need to sum to 1:
  Pj <- exp(logPj)
  Pj <- Pj/matrix(colSums(Pj), nrow = nrow(Pj), ncol = ncol(Pj), byrow = T)
  logPj <- log(Pj)

  return(logPj) # C by 2^K

}# end Neg_Binom_IRF()





#
Poisson_IRF <- function(param_j,  j, eta, X, categories.j_j){

  inter <- matrix(param_j[grep('inter', names(param_j))], nrow = 2, ncol = length(grep('inter', names(param_j))), byrow=T)
  slope <- matrix(c(0,1) * param_j['slope'], nrow = 2, ncol = ncol(inter), byrow = F)

  mu <- slope + inter
  mu <- exp(mu)

  mu <- matrix(mu, nrow = length(categories.j_j), ncol = nrow(mu), byrow = T) # Number of categories by number of latent classes
  Xj <- matrix(categories.j_j, nrow = length(categories.j_j), ncol = ncol(mu), byrow = F) # same, C by 2^K
  # Note that the Xj is NOT the item responses, it's just the possible categories
  # (NB: mu is usually written as lambda in textbooks)

  logPj <-  Xj * log(mu) - mu - lgamma(Xj + 1) # Note that this is NOT reduced to just the mean
  # This is C by 2^K

  # The probabilities need to sum to 1:
  Pj <- exp(logPj)
  Pj <- Pj/matrix(colSums(Pj), nrow = nrow(Pj), ncol = ncol(Pj), byrow = T)
  logPj <- log(Pj)

  return(logPj) # C by 2^K

}# end Poisson_IRF()



# Negative Binomial with Zero-Inflation Process
ZINB_IRF <- function(param_j,  j, eta, X, categories.j_j){

  inter <- matrix(param_j[grep('inter', names(param_j))], nrow = 2, ncol = length(grep('inter', names(param_j))), byrow=T)
  slope <- matrix(c(0,1) * param_j['slope'], nrow = 2, ncol = ncol(inter), byrow = F)

  mu <- slope + inter
  mu <- exp(mu) # ensures this parameter is positive
  size <- exp(param_j['size'])
  # 7.22.19, This is an easy way to constrain the size parameter
  # to be real & positive (see Bolker textbook pg 167 Ecological Models and Data in R)
  # Note that in the future you should then exponentiate the param_j['size'] because
  # that will be what the model actually used

  prob <- size/(size + mu)
  prob <- matrix(prob, nrow = length(categories.j_j), ncol = nrow(prob), byrow = T) # Number of categories by number of latent classes
  Xj <- matrix(categories.j_j, nrow = length(categories.j_j), ncol = ncol(prob), byrow = F) # same, C by 2^K
  # Note that the Xj is NOT the item responses, it's just the possible categories

  # Zero-inflation ADDED:
  infprob_j <- 1/(1 + exp(-1*(param_j['zi']))) #scalar
  # If response = 0
  #Can be 0 two different ways, either inflation OR from the count process
    # 'OR' means you add the probabilities (then take the log)
  prob0 <- log(infprob_j + (1-infprob_j)*prob^size)
  # If response > 0
  # multiply probability that it's NOT inflated zero times count process
  # This is an "AND" statement, requires multiplying probabilities
  prob1 <- log(1-infprob_j) +  lgamma(Xj + size) - lgamma(size) - lgamma(Xj + 1) + size*log(prob) + (Xj)*log(1-prob)

  logPj <- (Xj == 0) * prob0  +   (Xj != 0) *  prob1

  # The probabilities need to sum to 1:
  Pj <- exp(logPj)
  Pj <- Pj/matrix(colSums(Pj), nrow = nrow(Pj), ncol = ncol(Pj), byrow = T)
  logPj <- log(Pj)

  return(logPj) # C by 2^K

}# end ZINB_IRF()




# Zero Inflation Process
ZIP_IRF <- function(param_j,  j, eta, X, categories.j_j){

  inter <- matrix(param_j[grep('inter', names(param_j))], nrow = 2, ncol = length(grep('inter', names(param_j))), byrow=T)
  slope <- matrix(c(0,1) * param_j['slope'], nrow = 2, ncol = ncol(inter), byrow = F)

  mu <- slope + inter
  mu <- exp(mu) # ensures this parameter is positive
  mu <- matrix(mu, nrow = length(categories.j_j), ncol = nrow(mu), byrow = T) # Number of categories by number of latent classes
  Xj <- matrix(categories.j_j, nrow = length(categories.j_j), ncol = ncol(mu), byrow = F) # same, C by 2^K
  # Zero-inflation Probability:
  infprob_j <- 1/(1 + exp(-1*(param_j['zi']))) #scalar
  # If response = 0:
    # Hurdle model
    # logistic model for P(X=0)
    # Count model for P(X>0)
    # Zero-inflation models are NOT hurdle models
    # inflation model can include zeros, but they're distinct from the "excess" zeros
    # Two ways you can get zero, first term is "excess zero", second term is zero in Poisson process

  # If response = 0
    # Can be zero in one of two ways, either from the zero-inflation process OR from the count process
  prob0 <- log(infprob_j + (1-infprob_j)*exp(-mu))
  # If response > 0
    # Can be greater than zero IF you are not in the zero-inflation process AND you're in the count process
    # Note that 'AND' require multiplying probabilities (adding logs)
 # prob1 <- log(1-infprob_j) +  Xj * log(lambda_j) - lambda_j - lgamma(Xj + 1)
  prob1 <- log(1-infprob_j) +  Xj * log(mu) - mu - lgamma(Xj + 1)

  logPj <- (Xj == 0) * prob0  +   (Xj != 0) *  prob1


  # The probabilities need to sum to 1:
  Pj <- exp(logPj)
  Pj <- Pj/matrix(colSums(Pj), nrow = nrow(Pj), ncol = ncol(Pj), byrow = T)
  logPj <- log(Pj)

  return(logPj) # C by 2^K

}# end ZIP_IRF()



# DINA_IRF <- function(param_j, Xj, j, eta, N, K, alpha, Z = NULL, item.covar = NULL){
#
#   # specify the covariates to be included in item j:
#   if(is.null(item.covar[[j]])){Zj = NULL}else{Zj <- Z[ , item.covar[[j]], drop = F]}
#
#
#   p <- matrix(eta[ , j] * (1-param_j['slip']), N, 2^K, byrow = T) +
#           matrix((1-eta[ , j]) * param_j['guess'], N, 2^K, byrow = T)
#
#   llj <- Xj*log(p) + (1-Xj)*log(1-p)
#
#   return(llj)
#
# }# end DINA_IRF()


