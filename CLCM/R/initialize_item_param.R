
initialize_item_param <- function(item.type, X){


  J = length(item.type)
  param <- vector(mode="list", length=J)

for(j in 1:J){


#
  if(item.type[j] == 'Beta'){

    shape1 <- 3
    shape2 <- (shape1 - mean(X[ , j], na.rm = T))/mean(X[ , j], na.rm = T)
      #shape1 <- exp(shape1) # they must be positive! TODO: implement this constraint
    intercept <- 1
    param[[j]] <- c('slope' = shape1 - intercept, 'intercept' = intercept, 'shape2' = shape2)

  }# end Beta


#
  if(item.type[j] == 'Normal'){

    param[[j]] <- c(2 , mean(X[ , j], na.rm = T), sd(X[ , j], na.rm = T))
    names(param[[j]]) <- c('slope', 'intercept', 'sigma')

  }# end Normal


#
  if(item.type[j] == 'Nominal'){

    # p <- prop.table(table(X[ , j]))
    # xb <- -1*log( (1/p) - 1 )
    # all(p == exp(xb)/(1 + exp(xb))) # CHECK
    tmp <- rbind(rep(0, length(unique(X[ , j])) - 1),
                 rep(0.5, length(unique(X[ , j])) - 1)  )
    rownames(tmp) <- c('intercept', 'slope')
    colnames(tmp) <- paste0('k_', 1:ncol(tmp))
    param[[j]] <- tmp
    ## names(param[[j]]) <-

  }#end Ordinal


#
  if(item.type[j] == 'Ordinal'){

          tmp <- vector()
       for(score in (min(X[, j], na.rm =T) + 1):max(X[ ,j], na.rm =T)){ # 2.20.21 - start at 0
         tmp <- c(tmp, log(mean(X[ ,j] >= score, na.rm =T)))
       }

        param[[j]] <- c(2, tmp)
        names(param[[j]]) <- c('slope', paste0('intercept_', 1:(max(X[ ,j], na.rm =T))))

  }#end Ordinal



#
  if(item.type[j] == 'ZINB'){

    size = var(X[,j])/(mean(X[,j]) + mean(X[,j])^2)
    prob = size/(size + mean(X[,j]))
    zi = log(mean(X[ ,j] == 0, na.rm = T))
    param[[j]] <- c(2, prob, size, zi)
    names(param[[j]]) <- c('slope', 'intercept', 'size', 'zi')

  }#end ZINB



 #
  if(item.type[j] == 'ZIP'){

    param[[j]] <- c(2, -1*log(mean(X[ ,j][which(X[,j] != 0)])), log(mean(X[ ,j] == 0, na.rm = T)))
    names(param[[j]]) <- c('slope', 'intercept', 'zi')

  }#end ZIP



 #
  if(item.type[j] == 'Poisson'){

    param[[j]] <-c(2, log(mean(X[,j], na.rm = T)))
    names(param[[j]]) <- c('slope', 'intercept')

  }# end Poisson


 #
  if(item.type[j] == 'Neg_Binom'){

    size = var(X[,j])/(mean(X[,j]) + mean(X[,j])^2)
    prob = size/(size + mean(X[,j]))
    param[[j]] <- c(1, prob, size)
    names(param[[j]]) <- c('slope', 'intercept', 'size')

  }#end Neg Binom

}# end J


  # Item names
  if(!is.null(colnames(X))){
    names(param) <- colnames(X)
  } else {
    names(param) <- paste0('Item_', 1:J, '_', item.type)
  }


  # Initialize the number of item categories:
  categories.j <- vector(mode="list", length=J)

for(j in 1:J){

  categories.j[[j]] <-

    if(!(item.type[j] %in% c('Beta', 'Normal'))){

      min(X[, j], na.rm = T):max(X[, j], na.rm = T)

    } else {
      NA
    }
}# end loop over j


   return(list('param' = param, 'categories.j' = categories.j))


} #end initialization function

