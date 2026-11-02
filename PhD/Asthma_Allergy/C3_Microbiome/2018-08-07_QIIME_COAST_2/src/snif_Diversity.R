#
# similarity-sensitive diversity indices
# as per:
# slides: https://www.gla.ac.uk/media/media_349517_en.pdf
# slides: http://www.maths.ed.ac.uk/~tl/ncse/ncse_talk.pdf
# paper:  https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/10-2402.1
#         http://www.maths.ed.ac.uk/~tl/mdiss.pdf 
#
# NB (Alex): found bug in slide formula for q!=1: as per the paper p4, you actually need to take the reciproqual of the formula given.
# NB2 (Howard): discussed with Alex, formula in slide is correct. Not a bug. However may still give odd result. Howard to check code at some stage.

q.diversity <- function( abundance.mat, q=1, rel.mat=diag(rep(1), ncol(abundance.mat))  ) {

  # checks
  stopifnot( ncol(rel.mat)==ncol(abundance.mat) )
  stopifnot( nrow(rel.mat)==ncol(abundance.mat) )
  stopifnot( q>=0 )
  stopifnot( inherits(abundance.mat, "data.frame") | inherits(abundance.mat, "matrix")  )
  abundance.mat <- as.matrix(abundance.mat)
  stopifnot(isSymmetric.matrix(rel.mat)) # for now not supporting non-symmetric relationships 
  
  
  # get relative abundances
  prop.mat <- sweep( abundance.mat, MARGIN=1, FUN="*", STATS=1/rowSums(abundance.mat)  )
  
  
  # get ordinarinesses
  ord.mat <- prop.mat %*% rel.mat
  
  
  # calculate
  Dq <- rep(NA, nrow(prop.mat)) ; names(Dq) <- rownames(prop.mat)
  if (q==1) { 
    if (any(ord.mat==0)) stop( "Calculation will fail if there are any 0 ordinarinesses."  ) # TODO: drop species where ord is 0
    Dq <- exp(      -1 * apply(prop.mat * log( ord.mat ),           1, sum) )
  } else { 
    # for q!=1, use a generalised power mean
    Dq <- exp( log( 1 / apply(prop.mat *    ( ord.mat ^ (q - 1) ), 1, sum) ) / (q - 1) )
    # Dq <- rep(NA, nrow(prop.mat))
    # for (ind in 1:nrow(prop.mat)) {
    #   Dq[ind] <- 1 / (sum( prop.mat[ind,] * ((ord.mat[ind,])^(q-1)) ) ^ (1/(q-1)))
    # }
  }
  
  return(Dq)
}

q.diversity.profile <- function( abundance.mat, qs = seq(from=0, to=8, by=0.2), leg.pos="topright", ... ) {
  Dqs <- list()
  for (q in qs) {
    Dqs[[paste0("q.",q)]] <- q.diversity(abundance.mat, q, ...)
  }
  Dqs <- do.call(cbind, Dqs)
  if (is.null(rownames(Dqs))) rownames(Dqs) <- paste0(1:nrow(Dqs))
  plot( range(qs), range(Dqs), type="n", xlab="q", ylab="q-order diversity" )
  col.pal <- rainbow(nrow(Dqs))
  lapply(1:nrow(Dqs), function(i) {
    lines(  qs, Dqs[i,], col=col.pal[i] )
    points( qs, Dqs[i,], 
            col = ifelse(qs==1, col.pal[i], NA), 
            bg  = ifelse(qs==1, NA, col.pal[i]), 
            cex = ifelse(qs==1, 1, 0.5),
            pch = 21 )
  } )
  legend(leg.pos, legend = rownames(Dqs), text.col=col.pal)
  return(Dqs)
}

if (FALSE) {
  ## TEST ZONE
  abundance.mat <- matrix(c(1,24,6,1.5, 
                            1,12,3,2, 
                            2,19,6,2, 
                            0,13,5,1, 
                            1,14,6,2, 
                            5,5,5,5), byrow = TRUE, nrow=6, ncol=4) 
  colnames(abundance.mat) <- paste0("sp.",1:ncol(abundance.mat))
  rownames(abundance.mat) <- paste0("com.",1:nrow(abundance.mat))
  
  rel.mat <- diag(rep(1, ncol(abundance.mat))) ; rel.mat[2,1] <- rel.mat[1,2] <- 0.5 ; rel.mat[4,3] <- rel.mat[3,4] <- 0.25 ; rownames(rel.mat) <- colnames(rel.mat) <- colnames(abundance.mat)
  
  # q.diversity(abundance.mat, 0, rel.mat)
  # q.diversity(abundance.mat, 1, rel.mat)
  # q.diversity(abundance.mat, 2, rel.mat)
  
  q.diversity.profile(abundance.mat, rel.mat=rel.mat)
}


# EOF