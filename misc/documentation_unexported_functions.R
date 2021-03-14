

# Documentation for non-exported functions - pulled and deleted

# wtd.table()
#' Compute the weighted table for the Beta IRF
#' @param x pass the observed item responses to compute the ecdf
#' @param weights pass the posterior distribution to serve as weights
#' @return weighted ecdf values
#'
#'
#'
#' #' Likelihood
#'
#' Compute the loglikelihood of a CLCM for the computation of the jacobian. This computes
#' P(X) - probability of observing that item response pattern given the parameter estimates
#' This is used in the C2 function
#' Pass the same objects as "compute_loglike" PLUS add in the posterior

#' @param post posterior distribution, used to marginalize over the loglikelihood
#' @param param item parameters
#' @param X marix of item responses
#' @param item.type character vector of item types
#' @param eta matrix operationalizing the condensation rules
#' @param categories.j numeric vector of number of categories per item
#' @return matrix of N by 2^K


#' Prior Distribution
#'
#' Compute the log of the prior distribution for use in the EM algorithm
#'
#' @param post posterior distribution
#' @param Z data frame or matrix containing covariates
#' @param type the type of prior distribution - this is selected automatically in the clcm() function
#' depending on what is passed
#' @param reg.formula if a latent regression is estimated, pass the model specification as a character here
#' @return matrix of N by 2^K


#' Objective Function
#'
#' Compute the value of the objective function using vP parameter vector
#' @param vP vector of parameter values
#' @param X matrix of item responses
#' @param item.type_j  character vector of item type
#' @param eta matrix operationalizing the condensation rules
#' @param categories.j numeric vector of number of categories per item
#' @param support only used for Beta distribution item type, numeric vector that is the "support'
#' of the weighted ecdf used in the estimation routine.
#' @return matrix of N by 2^K


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


#' Loglikelihood
#'
#' Compute the loglikelihood of a CLCM
#' @param X matrix of item responses
#' @param item.type character vector of item types
#' @param eta matrix operationalizing the condensation rules
#' @param categories.j numeric vector of number of categories per item
#' @return matrix of N by 2^K

#' Item Parameters
#'
#' Initialize the item parameters for use in the CLCM estimation routine.
#' @param item.type character vector of item types
#' @param X matrix of item responses
#' @return List of item parameters


#' Item Probabilities
#'
#' Compute the item probabilities, latent class by category, for use in the E-step of the
#' EM algorithm
#'
#' @param item.type_j the type of item j
#' @param param_j the item parameters for item j
#' @param j item j
#' @param eta is the condensation rule
#' @param categories.j_j the number of categories in item j
#' @return Item probabilities, latent class by category


#' Alpha
#'
#' Create the Alpha pattern, the unique latent classes.
#' @param K the number of attributes or factors
#' @return 2^K by K matrix of 0/1 indicating the different latent classes


#' C2 Stage 1
#'
#' Compute the contrast matrices for the C2 Function
#' @param mod pass the clcm model
#' @param patt unique item response patterns
#' @param obs.cat.j unique observed categorical answers for item j
#' @return Contrast matrices to compute the C2 function


#' Pseudo Counts
#'
#' Computes the pseudo-counts or expected values as part of the E-step in the EM algorithm
#'
#'
#' @param post posterior distribution
#' @param X the item responses
#' @param j item number
#' @param categories.j_j the number of categories in the item (if applicable)
#' @param item.type_j the type of item
#' @return the expected number of subjects in each latent class responding in each category
