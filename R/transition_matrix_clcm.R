#' Transition Matrix
#'
#' Compute Transition Matrix, stratified on categorical variable
#' Key to evaluating group differences
#'
#' @param mod estimated model object from clcm() function.
#' Note that if estimating the transition matrix stratified on a covariate, that (categorical) covariate
#' must be part of the dataframe (dat) that was used to estimate the model,
#' i.e., mod$dat must contain the covariate
#'
#' @param modal.classification logical; classify subjects using modal a priori (MAP) classification?
#'
#' @param eap.classification select if expected a priori (EAP) classification is desired
#' If neither MAP nor EAP classification is selected, then sample-level averages will
#' be computed for each latent class. That is, the probablistic classifications in the
#' subject posterior distributions will be retained and averaged.
#'
#' @param threshold numeric value, if EAP classification is selected, must choose a threshold
#' for classification as 1 versus 0 on each attribute (factor).
#'
#' @param stratification logical; should the transition matrix be computed stratified on categorical covariate?
#'
#' @param covariate categorical variable, separate transition matrix estimated for each level of the variable.
#' Note that if estimating the transition matrix stratified on a covariate, that (categorical) covariate
#' must be part of the dataframe (dat) that was used to estimate the model,
#' i.e., mod$dat must contain the covariate.
#' @return Returns a 2^K by 2^k numeric matrix; if transition matrix is stratified on covariate,
#' then returns a list of 2^K by 2^K numeric matrices.
#' @export
#' @examples
#' \dontrun{
#'
#' set.seed(3112021)
#' sim.dat <- simulate_clcm(N=200,
#'                           number.timepoints = 2,
#'                           item.type = rep('Ordinal', 5),
#'                           categories.j = rep(4, 5),
#'                           lc.prop = list('Time_1' = c(0.5, 0.5), 'Time_2' = c(0.5, 0.5)) )
#'
#' mod <- clcm(dat = sim.dat$dat,
#'           item.type = sim.dat$item.type,
#'            item.names = sim.dat$item.names,
#'            Q = sim.dat$Q)
#'
#'
#' tau.hat <- transition_matrix_clcm(mod)
#'
#'
#' }
transition_matrix_clcm <- function(mod,
                                       eap.classification = F,
                                       threshold = NULL,
                                       modal.classification = F,
                                       stratification = F,
                                       #lat.reg = NULL){
                                       #Z = NULL,
                                       covariate = NULL){

 if(stratification == T & (is.null(covariate))){
    stop('Specify a covariate for stratification, for example: covariate = "Group")') }

 # if(eap.classification == T & is.null(threshold)){
  #  stop('If you want to use EAP classification, pass a threshold to use as a cut-off') }


  dat <- mod$dat
  alpha <- mod$alpha
  post.names <- mod$post.names
  K <- mod$K


  if(length(unique(dat$Time)) < 2) {
    stop('Cannot estimate transition matrix with fewer than two timepoints') }


  # if modal classification is used, overwrite the posteriors using the MAP posteriors:
  if(modal.classification == T){

      tp <- 1
    for(tt in unique(dat$Time) ){ # loop over timepoints

      ii <- which(dat[ , 'Time'] == tt)
      map <- matrix(0, nrow = length(ii), ncol = 2^K)
      map[cbind(1:N, apply(dat[ ii , post.names ], 1, which.max))] <- 1 # TODO: check the which.max tie breaker
      assign(x = paste0('post', tp), value = map)
      assign(x = paste0('Z', tp), value = dat[ ii, covariate, drop = F])
      tp <- tp + 1

      }

  } else { # end modal classification TRUE


      if(eap.classification == T){

          tp <- 1
        for(tt in unique(dat$Time) ){ # loop over timepoints

          ii <- which(dat[ , 'Time'] == tt)
          eap <- matrix(0, nrow = length(ii), ncol = 2^K)
          post.ii <- as.matrix(dat[ ii , post.names ])
          eap[cbind(1:N,  1 + (post.ii %*% alpha > threshold))]  <- 1 # EAP - customize threshold
          assign(x = paste0('post', tp), value = eap)
          assign(x = paste0('Z', tp), value = dat[ ii, covariate, drop = F])
          tp <- tp + 1

        }

      } else{  # end EAP classification TRUE


            # Not Modal, Not EAP, then it's just, whatever this is called:
                tp <- 1
            for(tt in unique(dat$Time) ){ # loop over timepoints

              ii <- which(dat[ , 'Time'] == tt)
              #eap <- matrix(0, nrow = length(ii), ncol = 2^K)
              post.ii <- as.matrix(dat[ ii , post.names ])
              #eap[cbind(1:N,  1 + (post.ii %*% alpha > threshold))]  <- 1 # EAP - customize threshold
              assign(x = paste0('post', tp), value = post.ii)
              assign(x = paste0('Z', tp), value = dat[ ii, covariate, drop = F])
              tp <- tp + 1

            }
      }
} # end ifelse statement on LCA used in computation of post1 and post2



  if(stratification == T){

      groups <- sort(unique(dat[ , covariate]))
      tau.hat <- vector(mode="list", length=length(groups))
      names(tau.hat) <- groups


      for(gg in groups){

          th <- (t(post1[which(Z1[ , covariate] == gg), ] ) %*% post2[which(Z2[ , covariate] == gg), ] ) /
                         matrix(colSums(post1[which(Z1[ , covariate] == gg), ]), 2^K, 2^K, byrow = F)


                  tau.hat[[paste0(gg)]] <- th

          }#end loop over covariates




    } else {


            tau.hat <- (t(post1)  %*% post2) /  matrix(colSums(post1), 2^K, 2^K, byrow = F)

        }



  return(tau.hat)

  } #END FUNCTION
