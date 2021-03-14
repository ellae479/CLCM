
C2_stage1 <- function(mod, patt, obs.cat.j){ # likely faster if only passing data.full[ , item.names]

  # Pull the following from the mod passed:
  dat.full <- mod$dat
  item.names <- mod$item.names
  dat <- dat.full[ , item.names]
    dat <- dat[complete.cases(dat), ] # CANNOT BE PERFORMED ON MISSING DATA
    dat <- lapply(dat[, colnames(dat)], factor)
    dat <- as.data.frame(dat)
  Y <- model.matrix( ~ . , data = dat)
    Y <- Y[ , -1]
    # N by (J times Kj) - dummy coded categories


# Unique Item Response Patterns:
  # mat <- apply(dat, 2, function(x) sort(unique(x)) )  # TODO does this need to be a factor or numeric?
  # mat <- expand.grid(mat, stringsAsFactors = T)
  # mat <- lapply(mat[, colnames(mat)], factor)
  # mat <- as.data.frame(mat)
  # mat <- model.matrix( ~ . , data = mat)
  # mat <- mat[ , -1] # drop intercept
  mat <- lapply(patt[, colnames(patt)], factor)
  mat <- as.data.frame(mat)
  mat <- model.matrix( ~ . , data = mat)
    mat <- mat[ , -1] # drop intercept

# Compute the proportion of subjects with each pattern:
  # I.E., Observed Probability for each unique item response pattern
    p.obs <- vector()
  for(i in 1:nrow(mat)){

    tmp <- ( Y == matrix(mat[i, ], nrow = nrow(Y), ncol = ncol(mat), byrow = T) )
    tmp <- apply(tmp, 1, function(x) sum(x) == ncol(Y))
    tmp <- mean(tmp)
    p.obs <- c(p.obs, tmp)

    }

  p.obs <- matrix(p.obs)
  # CHECK
  # sum(p.obs) # sums to 1, yeah this seems correct


# Create first order marginals
  L1 <- t(mat)
  p.obs1 <- L1  %*% p.obs  # First order marginals
  # CHECK
  # cbind(p.obs1, colMeans(Y)) # aligned

# Create the second order marginals:
# 1. Need rownames for the second order marginals - label L2 matrix to help clarify what you're looking at later:
  tmp <- mat[apply(mat, 1, sum) == 2, ]
  labels <- apply(t(apply(tmp, 1, function(x) colnames(tmp)[which(x == 1)])), 1, paste0, collapse = '')


    L2 <- vector()
    #loop over item pairs, i, using "which(apply(mat, 1, sum) == 2)"
    for(i in which(apply(mat, 1, sum) == 2)){

      # next, pull the unique response patterns where that item pair is 1, 1
      tmp1 <- mat * matrix(mat[i, ], nrow = nrow(mat), ncol = ncol(mat), byrow = T)
      # tmp 2 is which of the 128 unique response patterns has that item pair is 1, 1
      tmp2 <- apply(tmp1, 1, function(x) sum(x) == 2)
      L2 <- rbind(L2, tmp2)

     }


    L2 <- 1*L2
    rownames(L2) <- labels
    # dim(L2) -- r by c
    # "r" out of "c" response patterns have a unique pair of items that are 1,1

    # CHECK
    # p2 <- L2 %*% p.obs # second order marginals
    # cbind(p2, mat[apply(mat, 1, sum) == 2, ])
    # (t(Y) %*% Y)/nrow(Y)


#########################################
# Kronecker products

  dat <- dat.full[ , item.names]# Do this again so it's not a factor
  # CANNOT BE PERFORMED ON MISSING DATA
  dat <- dat[complete.cases(dat), ]
  all.comb <- combn(ncol(dat), 2)
  out <- list()


for(comb in 1:ncol(all.comb)){

  i <- all.comb[1, comb]
  j <- all.comb[2, comb]
  #tmp <- apply(dat, 2, function(x) sort(unique(x))[-1] )
  tmp <- lapply(obs.cat.j, function(x) x[-1]) # 3.6.21  Try this
  out.kr <- kronecker(X = tmp[[i]], Y = tmp[[j]])
  out[[paste0(i, j)]] <- out.kr

}


  out <- lapply(out, function(x) matrix(x, nrow = 1))
  R <- Matrix:::bdiag(out) # Create a block diagonal
  L.star <- R %*% L2
  row.names(L.star) <- apply(combn(colnames(dat), 2), 2, function(x) paste0(x[1], '_', x[2]))
  M <- rbind(L1, L.star)

#R %*% pi2
# 10x30 by 30x1
# "is a vector containing the I(I-1)/2 model-implied second order moments"
# pg 14 of the C2 manuscript
# TODO: check if this is correct, can it be greater than 1?
# Do the rows of pi2 correspond to the pairs output from combn function?
# Eyeballing it looks correct - confirm later


# Let M be a q x k matrix that vertically concatenates L1 and L.star such that
# its q1 first rows come from L1 and the remaining I(I-1)/2 rows come from L.star
# K is the number of observed probabilities (K = 128 here)
# I is the number of items (I = 5)
# D is the number of parameters


# pg 16 of the C2 manuscript
#M %*% pi.hat

out <- list('M' = M,
            'p.obs' = p.obs,
            'L1' = L1,
            'L2' = L2,
            'R' = R,
            'L.star' = L.star
            )

return(out)


}# end function
