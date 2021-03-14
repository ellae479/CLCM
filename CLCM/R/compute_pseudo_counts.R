
 pseudo_counts <- function(post, X, j, categories.j_j, eta, item.type_j){

    # library(Hmisc)
    post <- as.matrix(post)
    X <- as.matrix(X)

  if(item.type_j == 'Normal'){

     # Normal
    Xjc <- X[ , j, drop = F]
    post.reduced <- post %*% cbind(1- eta[ , j], eta[, j])
    c.l <- (t(Xjc) %*% post.reduced)/colSums(post.reduced, na.rm = T) # mean, stratified on LC
    SSE <- (matrix(Xjc, nrow = nrow(Xjc), ncol = ncol(post.reduced), byrow = F) -
                matrix(c.l, nrow = nrow(post.reduced), ncol = ncol(post.reduced), byrow = T) )^2
    SSE <- colSums(SSE * post.reduced, na.rm = T)
    # Pool variance
    deg.free <- colSums(post.reduced, na.rm = T) - 1  # n-1
    sigma <- sqrt(sum(SSE)/sum(deg.free))  # TODO resolve whether the deg free is appropriate: is it n-1?
    c.l <- c(c.l, sigma)
    # THis is the mean of latent class 1, mean of latent class 2, pooled standard deviation (pooled across latent classes)
    names(c.l) <- c(paste0('Mean_', c(0,1) ), 'sigma')


  } else {

    if(item.type_j == 'Beta'){


        # Beta
        post.reduced <- post %*% cbind(1- eta[ , j], eta[, j])
        support <- X[ , j, drop = F]
        c.l <- matrix(NA, nrow = length(support), ncol = ncol(post.reduced))

        for(l in 1:ncol(post.reduced)){
          ## w <- Hmisc::wtd.table(x = support, weights = post.reduced[ , l])
          w <- wtd.table(x = support, weights = post.reduced[ , l])
          cumu <- cumsum(w$sum.of.weights)
          c.l[ , l] <- cumu/cumu[length(cumu)]
        }



    }  else  {

            post.reduced <- post %*% cbind(1- eta[ , j], eta[, j])
            c.l <- matrix(NA, nrow = ncol(post.reduced), ncol = length(categories.j_j))

            for(cc in categories.j_j ){

              Xjc <- 1*(X[ , j, drop = F] == cc)
              c.l[ , cc + 1] <- t(Xjc) %*% post.reduced

            }# end loop over categories

            c.l <- t(c.l)  # Output should be: C by 2^K

  }} # end if else statement

  return(c.l)

} #end pseudo_counts

