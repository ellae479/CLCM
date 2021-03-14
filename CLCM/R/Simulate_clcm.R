#' Simulate Data
#'
#' Simulate data for the CLCM
#' @param N integer specifying the sample size
#' @param number.timepoints integer specify the number of timepoints, 1 or 2
#' @param Q the Q-matrix, a matrix of 1s and 0s specifying the factor loading structure.
#' The default is 1 factor (K=1), which forms two latent classes
#' @param item.type character vector specifying the type of item to be modeled
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
#'
#' @param categories.j numeric vector specifying the number of categories of each item.
#' For 'Normal' or 'Beta' item types, the value should be NA
#'
#' @param lc.prop list of the latent class proportions at each timepoint.
#' For example, `lc.prop = list(c('Time_1' = c(0, 1), 'Time_2' = c(0.6, 0.4))`
#' specifies that at timepoint 1 all subjects are in latent class 2,
#' and at timepoint 2, 40% are in latent class 2, with the remainder in latent class 1

#' @param transition.matrix a 2^K by 2^K numeric matrix that specifies the transition probabilities.
#' This is used in conjunction with the lc.prop at timepoint 1.
#' See Vignettes for a detailed example.
#'
#' @param post a matrix of the true posterior distributions - long data format.
#' Generate the posterior distributions according to user-preference, then pass posterior distributions
#' to the function.
#' This is the preferred way to specify the latent class proportions and transition probabilities because it
#' offers maximum control and flexibility.
#' See Vignettes for detailed examples on generating posterior distributions.
#'
#' @param param list of item parameters, default is to use the values in the function
#' @param item.names character vector of item names
#' @return Returns simulated item responses and posterior distributions.
#' @export
#' @examples
#' \dontrun{
#' set.seed(3112021)
#' simulate_clcm(N=50, number.timepoints = 1,
#'                item.type = rep('Ordinal', 5),
#'                categories.j = rep(4, 5),
#'                lc.prop = list('Time_1' = c(0.5, 0.5)) )
#' }
#'

simulate_clcm <- function(N,
                          number.timepoints,
                          Q = NULL,
                          item.type = NULL,
                          categories.j = NULL,
                          transition.matrix = NULL,
                          lc.prop = NULL,
                          post = NULL,
                          param = NULL,
                          item.names = NULL){

  # Check:
    if(is.null(item.type)){ stop('Specify Item Type - see documentation for options') }
    if(is.null(categories.j)){ stop('Specify the number of categories per item - see documentation for details') }
    if( is.null(lc.prop) & is.null(post) ){ stop('Pass latent class structure - see documentation for details') }
    if( !is.null(lc.prop) & !is.null(post) ){ stop('Pass **either** latent class proportions or matrices')}
    if(is.null(Q)){  print('No Q-matrix passed to estimation function; default is two latent classes');
                      Q <- matrix(1, nrow = length(item.type), ncol = 1, dimnames = list(paste0('Item_', 1:length(item.type)), NULL))
                      }

    #  Requires Q-matrix - create alpha and eta from Q-matrix
    # TODO extend eta to other types of condensation rules - currently only conjunctive
    # All condensation rules are equivalent for single attribute items
     K <- ncol(Q)
     alpha <- pattern(K)
     eta <- alpha %*% t(Q)
     eta <- ifelse(eta == matrix(1,2^K,1) %*% colSums(t(Q)),1, 0 )
     # pass eta to all functions


    dat <- data.frame(
                  'USUBJID' = rep(paste0('Subject_', formatC(1:N, width = 4, flag = '0')), length.out= N*number.timepoints),
                  'Time' = rep(paste0('Time_', 1:number.timepoints), each = N),
                  stringsAsFactors=F)


    # Initialize
    J <- length(item.type)

    # Likely will not pass item parameters, but is an option
    if(is.null(param)){

    param <- vector(mode = "list", length = J)
    names(param) <- paste0('Item_', 1:J, ' ', item.type)


  # Loop to set item parameter values:
    for(j in 1:J){


      if(item.type[j] == 'Ordinal'){
        tmp1 <- 2
            names(tmp1) <- 'slope'
        tmp2 <- seq(1, -2, length.out = categories.j[j] - 1)
            names(tmp2) <- paste0('intercept_', 1:(categories.j[j] - 1))
        tmp <- c(tmp1, tmp2)
        param[[j]] <- tmp
        }

      if(item.type[j] == 'Neg_Binom'){
        tmp <- c(2, 1, 1)
          names(tmp) <- c('slope', 'intercept', 'size')
        param[[j]] <- tmp
        }

      if(item.type[j] == 'Poisson'){
        tmp <- c(2, 1)
          names(tmp) <- c('slope', 'intercept')
        param[[j]] <- tmp
        }

       if(item.type[j] == 'ZINB'){
        tmp <- c(2, 1, 1, 1)
          names(tmp) <- c('slope', 'intercept', 'size', 'zi')
        param[[j]] <- tmp
        }

      if(item.type[j] == 'ZIP'){
        tmp <- c(2, 1, 1)
          names(tmp) <- c('slope', 'intercept', 'zi')
        param[[j]] <- tmp
      }

      if(item.type[j] == 'Normal'){
        tmp <- c(2, 0, 1)
          names(tmp) <- c('slope', 'intercept', 'sigma')
        param[[j]] <- tmp
      }

      if(item.type[j] == 'Beta'){
        tmp <- c(2, 2, 2)
          names(tmp) <- c('slope', 'intercept', 'shape2')
        param[[j]] <- tmp
      }

      if(item.type[j] == 'Nominal'){
        # Matrix
        tmp1 <- seq(-2, -1, length.out = categories.j[j] - 1)
        tmp2 <- seq(1, 2, length.out = categories.j[j] - 1)
        tmp <- rbind(tmp1, tmp2)
        rownames(tmp) <- c('intercept', 'slope')
        colnames(tmp) <- paste0('k_', 1:(categories.j[j]-1))
        param[[j]] <- tmp

        }


    }# end loop to set parameter values


    } else {  # end if else for optional item parameter - likely will require generation in the above

      param <- param # if you pass item parameters,just use those

    }

    # Joint Distribution:
    # If passed the latent class proportions, then sample from those proportions

if(!is.null(post)){

  # if passed a matrix of subject-level latent classifications, use those
  # Note that can use latent regression to generate these - see Vignettes
  post <- post


  } else {

    if(!is.null(transition.matrix)){

          post <- matrix(0, nrow = nrow(dat), ncol = 2 ^ K)
              # First time point:
              tt <- sort(unique(dat$Time))[1]
              ii <- which(dat[ , 'Time'] == tt)
              post.1 <- post[ ii, ]
              post.1[ cbind(1:length(ii), sample( c(1:2^K), size = length(ii), replace = T, prob = lc.prop[[ tt ]]))] <- 1
              post[ ii ,  ] <- post.1
              # Second time point:
              tt <- sort(unique(dat$Time))[2]
              ii <- which(dat[ , 'Time'] == tt)
              post.2 <- post.1 %*% tau
              # Posterior Distributions all sum to 1?
                if(!all(apply(post.2, 1, sum) == 1)){ stop("Posterior distributions don't sum to 1; check latent class proportions and transition matrix") }
              p2 <- t(apply(post.2, 1, cumsum))
              lca2 <- apply(matrix(runif(n = nrow(p2)), nrow = nrow(p2), ncol = ncol(p2), byrow = T) > p2, 1, sum) + 1
              post.2[] <- 0
              post.2[ cbind(1:length(lca2), lca2) ] <- 1
              post[ ii ,  ] <- post.2

    } else {

              post <- matrix(0, nrow = nrow(dat), ncol = 2 ^ K)
          for(tt in unique(dat$Time) ){

              ii <- which(dat[ , 'Time'] == tt)
              post.t <- post[ ii, ]
              post.t[ cbind(1:length(ii), sample( c(1:2^K), size = length(ii), replace = T, prob = lc.prop[[ tt ]]))] <- 1
              post[ ii ,  ] <- post.t

          }# end tt

    }
  } # end Generation  of Joint Distribution



    # Initialize Item Responses:
    X <- matrix(NA, nrow = nrow(dat), ncol = J)

# Generate Item Responses:
  j <- 1
  for(j in 1:J){

    # Reduced Posterior Distribution
    # from 2^K down to 2, conjunctive models
    # use eta_j to create it
    post.reduced <- post %*% cbind(1- eta[ , j], eta[, j])

      if(item.type[j] == 'Beta'){
          Pj <- Return_logPj(item.type_j = item.type[j], param_j = param[[j]], j = j, eta = eta, categories.j_j = NA)
            shape2 <- Pj['shape2']
            Pj <- matrix(Pj[grep('shape1', names(Pj))], ncol = 1)
            ## shape1 <- post %*% Pj
            shape1 <- post.reduced %*% Pj
            X[ , j] <- rbeta(n = nrow(shape1), shape1 = shape1, shape2 = shape2)
      }# End Beta


      if(item.type[j] == 'Normal'){
          Pj <- Return_logPj(item.type_j = item.type[j], param_j = param[[j]], j = j, eta = eta, categories.j_j = NA)
            sigma <- Pj['sigma']
            Pj <- matrix(Pj[grep('mean', names(Pj))],  ncol = 1)
            ## mean.l <- post %*% Pj
            mean.l <- post.reduced %*% Pj
            X[, j] <- rnorm(n = nrow(mean.l), mean = mean.l, sd = sigma) #
      }# End Normal


      if(item.type[j] == 'Ordinal'){
          logPj <- Return_logPj(item.type_j = item.type[j], param_j = param[[j]], j = j, eta = eta, categories.j_j = NA)
            Pj <- exp(logPj)
            ## tmp <- post %*% t(Pj)
            tmp <- post.reduced %*% t(Pj)
            X[, j] <- apply(runif(n = N) > t(apply(tmp, 1, cumsum)), 1, sum) #+ 1  TODO: TEST! 2.20.21
      }# End Ordinal


      if(item.type[j] %in% c('Nominal', 'Neg_Binom', 'Poisson', 'ZIP', 'ZINB') ){
          logPj <- Return_logPj(item.type_j = item.type[j], param_j = param[[j]], j = j, eta = eta,
                                categories.j_j = 0:(categories.j[j]-1))
            Pj <- exp(logPj)
            ## Pij <- post %*% t(Pj)
            Pij <- post.reduced %*% t(Pj)
          for(i in 1:nrow(Pij)){
            X[i, j] <- sample(x = 0:(categories.j[j]-1), size = 1, prob = Pij[i, ])
          }
      }# End Neg Binom, Poisson, ZIP, ZINB


  }# End Item Generation Loop


  # Item Names
  if(is.null(item.names)){  item.names <- paste0('Item_', 1:J) }
    colnames(X) <- item.names
  post.names <- paste0('true_post_LC_', apply(alpha, 1, paste0, collapse = '')) # 3.6.21, fixed
    colnames(post) <- post.names

  lca <- apply(post, 1, which.max)
    #names(lca) <- 'true_lca'

  dat <- cbind.data.frame(dat, X, post, 'true_lca' = lca) # 3.8.21: added in lca for ease


    return(list('dat' = dat,
                'item.responses' = X,
                'post' = post,
                'lca' = lca,
                'param' = param,
                'item.type' = item.type,
                'item.names' = item.names,
                'categories' = categories.j,
                'lc.prop' = lc.prop,
                'Q' = Q,
                'alpha' = alpha,
                'K' = K,
                'eta' = eta))

} #end Simulate_Data function
