#' Fit CLCM
#'
#' Estimate confirmatory latent class model
#' @param dat data frame containing item responses. If the data contains more than
#' one timepoint, it must be specified by a variable named `Time`, a factor of two levels.
#' Data should be in long format, that is, similar format to that used by the nlme/lme4/glmmTMB R packages
#' @param Q optional pass the Q-matrix, default is K=1, two latent classes.
#' This is a confirmatory latent class model, hence the need for the specifcation.
#' Note that only dichotomous attributes (factors) are supported.
#' Only the conjunctive condensation rule is supported.
#' All condensation rules are equivalent for items that evaluate a single attriute (factor).
#' @param item.type character vector specifying the type of item to be modeled. The item type options
#' are as follows:
#' \itemize{
#' \item `Ordinal` - Ordinal model, 1 slope parameter, C-1 intercept parameters,
#' where C is the number of categories. The slope parameter distinguishes the latent classes.
#' \item `Nominal` - Nominal or Multinomial model, 2*(C-1) parameters. The slope parameters
#' distinguish the latent classes.
#' \item `Poisson` - Poisson model, 2 parameters.
#' The model is paramterized using lambda.
#' See ?rpois for details.
#' The lambda parameter is stratified on the latent classes.
#' The slope parameter distinguishes the latent classes.
#'
#' \item `Neg_Binom` - Negative Binomial model, 3 parameters.
#' The model is parameterized using mu and size.
#' See ?rnbinom for details on parameters.
#' The mu parameter is stratified on the latent classes; the size parameter is common to
#' all latent classes.
#'The slope parameter distinguishes the latent classes.
#'
#'
#' \item `ZINB` - Zero-Inflated Negative Binomial model, 4 parameters.
#' The model is parameterized as mu, size, and zi (zero inflation parameter).
#' The mu parameter is stratified on the latent classes;
#' the size and zi parameters are common to all latent classes.
#' The slope parameter distinguishes the latent classes.
#'
#'
#' \item `ZIP`- Zero-Inflated Poisson model, 3 parameters.
#' The model is parameterized using lambda and zi (zero inflation parameter).
#' The lambda parameter is stratified on latent classes.
#' The zi parameter is common to all latent classes.
#' The slope parameter distinguishes the latent classes.
#'
#'
#' \item `Normal` - Normal distribution model, 2 parameters.
#' The model is parameterized using a mean and variance
#' The mean is stratified on latent classes.
#' The variance is common to all latent classes (pooled across latent classes)
#' The slope parameter distinguishes the latent classes.
#'
#'
#' \item `Beta` - Beta distribution model, 3 parameters.
#' Support ranges from 0 to 1.
#' The model is parameterized using shape1 and shape2 paramters.
#' See ?rbeta for details
#' The shape1 parameter is stratified on latent classes
#' The shape 2 parameter is common to all latent classes
#' The slope parameter distinguishes the latent classes.
#' }
#' @param item.names specify the item names; the dataframe column names containing item responses
#'
#' @param lc.con optional list of constraints on each latent class at each timepoint, where
#' `NA` constrains that latent class to zero at that timepoint.
#' For example, constraint latent class one (out of two) to be zero at timepoint 1:
#' `lc.con = list('Time_1' = c(NA, 1), 'Time_2' = c(1, 1))`.
#' Note that the `dat` dataframe passed must contain the
#' variable `Time`.
#'
#'
#' @param lat.reg optional pass variables to regress latent class on
#' For example, regress the latent classes at timepoint 2 onto the Group variable:
#' `lat.reg = list('Time_1' = NULL, 'Time_2' = 'Group' )`, where `Time_1` and `Time_2` are
#' the discrete timepoints. Note that the `dat` dataframe passed must contain the
#' variable `Time`.
#'
#' @param sv optional list of the starting values for item parameter estimation
#' @param initial.post optional matrix of the initial posterior distribution
#' @param max.diff convergence tolerance of item param estimation; default is 1e-04
#' @param max.it maximum number of iterations in EM estimation procedure; default is 1e3
#' @param verbose logical; should the function print the estimation progress? Default is TRUE,
#' recommend set to FALSE for simulations.
#' @return Returns the estimated confirmatory latent class model
#' @export
#' @examples
#' \dontrun{
#' set.seed(3112021)
#' sim.dat <- simulate_clcm(N=200,
#'                           number.timepoints = 1,
#'                           item.type = rep('Ordinal', 5),
#'                           categories.j = rep(4, 5),
#'                           lc.prop = list('Time_1' = c(0.5, 0.5)) )
#'
#' mod <- clcm(dat = sim.dat$dat,
#'             item.type = sim.dat$item.type,
#'             item.names = sim.dat$item.names,
#'             Q = sim.dat$Q)
#'
#'  }


clcm <- function(dat,
                 item.type = NULL,
                 item.names = NULL,
                 Q = NULL,
                 lc.con = NULL,
                 lat.reg = NULL,
                 sv = NULL,
                 post.true = NULL,
                 initial.post = NULL,
                 max.diff = 1e-4,
                 max.it = 1e3,
                 verbose =T ){


# Check Data - Throw Errors, set defaults
  if(!is.data.frame(dat)){ stop('Pass a dataframe containing item responses')}
  if(is.null(item.type)){ stop('Must specify item type; pass character vector')}
  if(is.null(item.names)){ stop('Must specify item names; pass character vector')}
  if(is.null(dat$Time)){ print('No time variable - assuming a single timepoint');
                          dat$Time <- 'Time_1' }
  if(length(unique(dat$Time)) > 2){ stop('Estimation function currently only supports up to 2 timepoints; please check your "Time" variable to confirm you only have two timepoints')}
  if(is.null(Q)){  print('No Q-matrix passed to estimation function; default is two latent classes');
                      Q <- matrix(1, nrow = length(item.type), ncol = 1, dimnames = list(paste0('Item_', 1:length(item.type)), NULL))
  }

  # Set default values:
  # if(is.null(max.it)){  max.it <- 1000}
  # if(is.null(max.diff)){  max.diff <- 0.0001 }

  #  Requires Q-matrix - create alpha and eta from Q-matrix
  # TODO extend eta to other types of condensation rules - currently only conjunctive
  # All condensation rules are equivalent for single attribute items
     K <- ncol(Q)
     alpha <- pattern(K)
     eta <- alpha %*% t(Q)
     eta <- ifelse(eta == matrix(1,2^K,1) %*% colSums(t(Q)),1, 0 )
     # pass eta to all functions

  if(is.null(lc.con)){  lc.con <- list('Time_1' = rep(1, 2^K), 'Time_2' = rep(1, 2^K) )  }
  if(is.null(lat.reg)){  lat.reg <- list('Time_1' = NULL, 'Time_2' = NULL)  }

  #  Initialize Item Parameters
  ip <- initialize_item_param(item.type = item.type, X  = dat[ , item.names]  )
  categories.j <- ip[['categories.j']]
  if(!is.null(sv)){  param <- sv } else {  param <- ip[['param']]   }


# Initialize Priors - full data frame dat
  tmp <- matrix(NA, nrow = nrow(dat), ncol = 2^K)
  lprior.names <- paste0('lprior_LC_', apply(alpha, 1, paste0, collapse = '')) # 3.6.21, fixed
  colnames(tmp) <- lprior.names
  dat <- cbind.data.frame(dat, tmp)  # added columns for the prior here


  for(tt in unique(dat$Time) ){ # loop over timepoints

    lprior <- rep(1, 2^K)
    lprior[is.na( lc.con[[ tt ]] ) ] <- 0
    lprior <- lprior/sum(lprior)
    lprior <- log(lprior)

    ii <- which(dat[ , 'Time'] == tt)
    lprior <- matrix(lprior, nrow = length(ii) , ncol = length(lprior), byrow = T)
    dat[ ii , lprior.names ] <- lprior

    }


# Initialize Posterior Distribution - full data frame dat
  tmp <- matrix(NA, nrow = nrow(dat), ncol = 2^K)
  post.names <- paste0('post_LC_', apply(alpha, 1, paste0, collapse = '')) # 3.6.21, fixed
  colnames(tmp) <- post.names
  dat <- cbind.data.frame(dat, tmp)  # added columns for the prior here

# Initialize Posterior Distribution:
  post.full <- compute_post(X = dat[ , item.names],
                               lprior = dat[ , lprior.names], # Need labels from lprior initialization
                               item.type = item.type,
                               param = param,
                               eta= eta,
                               categories.j = categories.j)

  dat[ , post.names] <- post.full

  if(!is.null(initial.post)){ dat[ , post.names] <- initial.post }

####################################################################################################
#
#
#             START ESTIMATION
#
#
######################################################################################################3


  diff = 1
  it = 0

#Start estimation:
while(diff > max.diff & it < max.it){ #

  est0 <- param # use this to check item convergence later

 # Loop over J Items:
 for(j in 1:length(item.type) ){


  #E step:
  item.type_j <- item.type[[j]]
  categories.j_j <- categories.j[[j]]

  # Pass "post.names" and "item.names"
  ev = pseudo_counts(post = dat[ , post.names], X = dat[ , item.names],
                          j = j,  categories.j = categories.j_j,
                          eta = eta, item.type_j = item.type_j)


  # M step:
  vP <- param[[j]]

  suppressWarnings(

    out <- optim(par = vP, fn = objective_function, method = 'BFGS',
                 item.type_j = item.type_j,
                 eta = eta,
                 j = j,
                 ev = ev,
                 categories.j = categories.j_j,
                 support = dat[ , item.names])

  )
  # Update parameter estimates:
  param[[j]] <- out$par


}# end loop over J items


  # ########################
  # UPDATE POSTERIOR DISTRIBUTION
  post.full <- compute_post(X = dat[ , item.names],
                               lprior = dat[ , lprior.names], # Need labels from lprior initialization
                               item.type = item.type,
                               param = param,
                               eta= eta,
                               categories.j = categories.j)

  dat[ , post.names] <- post.full

  if(!is.null(post.true)){  dat[ , post.names] <- post.true  }


#############################################
  # UPDATE PRIOR DISTRIBUTION
for(tt in unique(dat$Time) ){ # loop over timepoints

  ii <- which(dat[ , 'Time'] == tt)
  post.t <- dat[ ii , post.names ]
  Z <- dat[ ii, ]

  if(!is.null( lat.reg[[ tt ]] )){

        lprior.t <- compute_lprior(post = post.t, K = K,  Z = Z, type = 'latent_regression', reg.formula = lat.reg[[ tt]] )
          lat.reg.param <- lprior.t$vP
          lprior.t <- lprior.t$lprior # outputs both the prior and the latent regression parameter estimates

     } else {

        lprior.t <- compute_lprior(post = post.t, K = K, type = 'standard_EB')

  }# end else statement

  dat[ ii , lprior.names ] <- lprior.t

}# end loop over timepoints




  # Convergence: item parameter estimates - note, can track the convergence history
  diff <- max(abs(unlist(est0) - unlist(param)))
  it = it + 1
  if(verbose == T){  cat('iteration: ', it, '  max diff in item parameter estimates: ', round(diff, 6), '\n') }

} # End Estimation "while" loop


  # Post Processing:
  lca <- apply(dat[ , post.names], 1, which.max)
  dat <- cbind.data.frame(dat, lca)

  if(!exists('lat.reg.param')) { lat.reg.param <- NULL }


  # Returns:
  out <- list('item.param' = param,
              'lat.reg.param' = lat.reg.param,
              'dat' = dat,
              'item.type' = item.type,
              'categories.j' = categories.j,
              'item.names' = item.names,
              'lprior.names' = lprior.names,
              'post.names' = post.names,
              'eta' = eta,
              'alpha' = alpha,
              'K' = K,
              'Q' = Q,
              'lc.con' = lc.con,
              'lat.reg' = lat.reg)

  return(out)

} #end function Estimate clcm

