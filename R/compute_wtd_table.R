
wtd.table <- function(x, weights=NULL, type=c('list','table'),
                      normwt=FALSE, na.rm=TRUE)
{

  # wtd.table() Function from the Hmisc package
  # https://github.com/harrelfe/Hmisc/blob/master/R/wtd.stats.s
  # Do not need the entire Hmisc R package, just this function

  type <- match.arg(type)
  if(! length(weights))
    weights <- rep(1, length(x))

  ## isdate <- testDateTime(x)  ## 31aug02 + next 2
  isdate <- F
  ax <- attributes(x)
  ax$names <- NULL

  if(is.character(x)) x <- as.factor(x)
  lev <- levels(x)
  x <- unclass(x)

  if(na.rm) {
    s <- ! is.na(x + weights)
    x <- x[s, drop=FALSE]    ## drop is for factor class
    weights <- weights[s]
  }

  n <- length(x)
  if(normwt)
    weights <- weights * length(x) / sum(weights)

  i <- order(x)  # R does not preserve levels here
  x <- x[i]; weights <- weights[i]

  if(anyDuplicated(x)) {  ## diff(x) == 0 faster but doesn't handle Inf
    weights <- tapply(weights, x, sum)
    if(length(lev)) {
      levused <- lev[sort(unique(x))]
      if((length(weights) > length(levused)) &&
         any(is.na(weights)))
        weights <- weights[! is.na(weights)]

      if(length(weights) != length(levused))
        stop('program logic error')

      names(weights) <- levused
    }

    if(! length(names(weights)))
      stop('program logic error')

    if(type=='table')
      return(weights)

    x <- all.is.numeric(names(weights), 'vector')
    if(isdate)
      attributes(x) <- c(attributes(x),ax)

    names(weights) <- NULL
    return(list(x=x, sum.of.weights=weights))
  }

  xx <- x
  if(isdate)
    attributes(xx) <- c(attributes(xx),ax)

  if(type=='list')
    list(x=if(length(lev))lev[x]
           else xx,
         sum.of.weights=weights)
  else {
    names(weights) <- if(length(lev)) lev[x]
                      else xx
    weights
  }
}
