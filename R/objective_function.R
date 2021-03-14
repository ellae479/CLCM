
objective_function <- function(vP, item.type_j, j, eta, ev, categories.j_j, support = NULL){


    logPj <- Return_logPj(param_j = vP, item.type_j = item.type_j, j=j, eta = eta, categories.j_j = categories.j_j)


    if(item.type_j == 'Normal'){

      LL <- sum( (ev - logPj)^2 , na.rm = F)  # sum squared diff between expected vs observed

    } else {

        if(item.type_j == 'Beta'){

          # Expected (Model Implied)
          support <- sort(support[ , j])
          shape1.l <- logPj[grep('shape1', names(logPj))] # shape1, for each latent class
          shape2 <- logPj['shape2'] # shape2 - not differentiated across latent classes - legit??
          e <- matrix(NA, nrow = length(support), ncol = length(shape1.l)) #Expected (model implied)
          for(i in 1:nrow(e)){
            e[i, ] <- pbeta(support[i], shape1 = shape1.l, shape2 = shape2, lower.tail = TRUE)
          } # end i loop

          LL <- sum( (ev - e)^2, na.rm = F)

        } else {

            LL <- ev * logPj
            LL <- -1*sum(LL, na.rm = F) #leaving in the NA is important

    }} #end if else statement


      return(LL)

}#end objective_function

