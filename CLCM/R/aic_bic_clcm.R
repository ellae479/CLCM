#' AIC, BIC, -2LL
#'
#' Compute AIC, BIC, and 2LL of the CLCM model
#' @param mod object from the `clcm()` function
#' @return Returns the -2LL, AIC, and BIC of the confirmatory latent class model.
#' @export
#' @examples
#' \dontrun{
#' set.seed(3112021)
#' sim.dat <- simulate_clcm(N=1000,
#'                          number.timepoints = 1,
#'                          item.type = rep('Ordinal', 5),
#'                          categories.j = rep(4, 5),
#'                          lc.prop = list('Time_1' = c(0.5, 0.5)) )
#'
#'
#' mod1 <- clcm(dat = sim.dat$dat,
#'             item.type = rep('Nominal', 5),
#'             item.names = sim.dat$item.names,
#'             Q = sim.dat$Q)
#'
#'
#'mod2 <- clcm(dat = sim.dat$dat,
#'            item.type = rep('Ordinal', 5),
#'            item.names = sim.dat$item.names,
#'            Q = sim.dat$Q)
#'
#'
#' mod.fit1 <- aic_bic_clcm(mod1)
#' mod.fit2 <- aic_bic_clcm(mod2)
#'
#'
#' unlist(mod.fit1)
#' unlist(mod.fit2)
#' }





aic_bic_clcm <- function(mod){  # pass the model output from "clcm" function

    dat <- mod$dat
    lprior.names <- mod$lprior.names
      lprior <- as.matrix(mod$dat[ , lprior.names])
    item.names <- mod$item.names
    item.type <- mod$item.type
    item.param <- mod$item.param
    eta <- mod$eta
    categories.j <- mod$categories.j
    lc.con <- mod$lc.con
    lat.reg <- mod$lat.reg
    lat.reg.param <- mod$lat.reg.param


    LL <- compute_loglike(X = dat[ , item.names],
                                  item.type = item.type,
                                  param = item.param,
                                  eta= eta,
                                  categories.j = categories.j)


    LL = LL + lprior
    LL = sum(log(rowSums(exp(LL), na.rm = T)), na.rm = T)

    # Number of parameters
      npar <- vector()
    for(tt in unique(dat$Time) ){

      if( !is.null(lat.reg[[tt]]) ) {

        npar <- c(npar, length(lat.reg.param) )

          } else {

            ## npar <- c(npar, length(lc.con[[tt]]) - 1)
            npar <- c(npar, sum(!is.na(lc.con[[tt]])) - 1) # 3.10.21: you have to pass a full list of constraints for this to work

          }
    }# end loop over timepoints

  npar <- sum(npar) # latent Class/Structural Model Parameters
  npar <- npar + length( unlist(item.param) )

  # TODO what is the "N' in bic with a longitudinal repeated measures model?
  N <- length(unique(dat$USUBJID))

  aic <- -2*LL + 2*npar
  bic <- -2*LL + log(N)*npar

  return(list('neg_2LL' = -2*LL, 'npar' = npar, 'AIC' = aic, 'BIC' = bic))

}#end compute_aic_bic()

#
