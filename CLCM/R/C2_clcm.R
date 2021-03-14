#' C2 Fit Statistic
#'
#' Compute C2 Test Statistic for the CLCM
#' @param mod estimated model object from `clcm()` function
#' @param verbose logical; print the C2 statistic and associated p-value? Default is TRUE, set to FALSE for simulations.
#' @return Returns the C2 fit statistic, associated p-value (test of exact fit,
#' where null hypothesis is the model fits the data exactly), RMSEA value,
#' sample size N, degrees of freedom, observed probabilities, and model-implied probabilities.
#' The latter two can be used to assess the fit visually.
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
#' mod.fit <- C2_clcm(mod)
#'
#'  }


C2_clcm <- function(mod, verbose = T){


  # Pull the following from the mod passed:
  dat.full <- mod$dat
  item.names <- mod$item.names
  dat <- dat.full[ , item.names]
  post <- dat.full[ , mod$post.names] # This is dumb, just pass the fucking posterior distributions dude
  eta <- mod$eta
  item.type <- mod$item.type
  categories.j <- mod$categories.j
  item.param <- mod$item.param

  if(!all(item.type == 'Ordinal')){ stop('The C2 statistic only works for ordinal items') }

# CANNOT BE PERFORMED ON MISSING DATA
  dat <- dat[complete.cases(dat), ]
# The unique item responses are the OBSERVED categories for each item
  # slightly diff from categories.j; this is possible categories
  # Ehh I guess for ordinal it's the same...shit
  obs.cat.j <- vector(mode="list", length=length(item.names))
  for(j in 1:length(item.names)){
    obs.cat.j[[j]] <- sort(unique(dat[ , j]))
  }
  patt <- expand.grid(obs.cat.j, stringsAsFactors = T)

  # 3.6.21: Replaced the below with the above
  # Create the unique item response patterns:
  # dat <- lapply(dat[, colnames(dat)], factor)
  # dat <- as.data.frame(dat)
  # Y <- apply(dat, 2, function(x) sort(unique(x)) )
  # Y <- expand.grid(Y, stringsAsFactors = T)
  # N <- nrow(Y)  # so N is not the fucking sample size???
  # #J <- ncol(Y) # Do we need this?



  loglike <- compute_loglike(X = patt, item.type = item.type, param = item.param, eta =eta, categories.j = categories.j)
  like <- exp(loglike)
  p.mod <-  like %*% colMeans(post) # Model Implied Probabilities

  # Compute Jacobian using numerical function:
      vP <- unlist(item.param)
    j.hat.full <- numDeriv::jacobian(func = compute_like_jacobian,
                           x = vP,
                           X = patt,
                           item.type = item.type,
                           eta = eta,
                           categories.j = categories.j,
                           post = post,
                           item.names = item.names)

  # Create the Reduction Matrices and Observed Probabilities
  stage1 <- C2_stage1(mod = mod, patt = patt, obs.cat.j = obs.cat.j)
  M <- stage1$M
    M <- as.matrix(M)
  j.hat <- M %*% j.hat.full
  q <- nrow(j.hat)
  d <- ncol(j.hat)
  # Jacobian hat, dimension q x d
  # q <- number of marginals being evaluated, equal to nrow(M)
  # d <- number of parameters estimated

  Xi <- diag(as.vector(p.mod)) - p.mod %*% t(p.mod)
# formula is from middle of page 15 of the C2 manuscript, see it there
  y.hat <- M %*% Xi %*% t(M)
  omega <- solve(y.hat) -
    solve(y.hat) %*% j.hat %*% (solve( t(j.hat) %*% solve(y.hat) %*% j.hat ) ) %*% t(j.hat) %*% solve(y.hat)

  #
  p.obs <- stage1$p.obs
  ## N <- nrow(dat) # 3.6.21: Type I error inflation, so switch to the N below
  N <- length(unique(dat.full$USUBJID)) # Option to change the N;
  df <- q - d
  C2 <- N * t(M %*% (p.obs - p.mod)) %*% omega %*% M %*% (p.obs - p.mod)


  # Final output:
  pval <- 1 - pchisq(C2, df = df)
  Fhat <- C2/N
  rmsea <- sqrt(Fhat/df)

  print.out <- c(C2, pval, rmsea)
  if(verbose == T){  cat(paste0(c('C2 = ','  P-Value = ', '  RMSEA = '),  sprintf("%0.3f", print.out) ), '\n')  }

  out <- list('C2' = C2,
              'pval' = pval,
             # 'Fhat' = Fhat,
              'rmsea' = rmsea,
              'N' = N,
              'df' = df,
              'p.obs' = p.obs,
              'p.mod' = p.mod )

    return(out)


}

