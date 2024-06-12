part_coeff <- function(g, memb, A=NULL, weighted=FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse=FALSE, names=FALSE, attr='weight')
    } else {
      A <- as_adj(g, sparse=FALSE, names=FALSE)
    }
  }
  Ki <- colSums(A)
  N <- max(memb)
  Kis <- t(rowsum(A, memb))
  1 - ((1 / Ki^2) * rowSums(Kis^2))
}

within_module_deg_z_score <- function(g, memb, A=NULL, weighted=FALSE) {
  stopifnot(is_igraph(g))
  if (is.null(A)) {
    if (isTRUE(weighted)) {
      A <- as_adj(g, sparse=FALSE, names=FALSE, attr='weight') #converts graph to adjacency matrix
    } else {
      A <- as_adj(g, sparse=FALSE, names=FALSE)
    }
  }
  N <- max(memb)
  nS <- tabulate(memb) #how many nodes are in each module
  z <- Ki <- rep.int(0, dim(A)[1L])
  Ksi <- sigKsi <- rep.int(0, N)
  
  for (S in seq_len(N)) {
    x <- rowSums(A[memb == S, memb == S, drop=FALSE])
    Ki[memb == S] <- x
    Ksi[S] <- sum(x) / nS[S]
    sigKsi[S] <- sqrt(sum((x - Ksi[S])^2) / (nS[S] - 1L))
  }
  z <- (Ki - Ksi[memb]) / sigKsi[memb]
  z[is.infinite(z)] <- 0
  return(z)
}

assign_module_roles <- function(zp){
  zp <- na.omit(zp)
  zp$roles <- rep(0, dim(zp)[1])
  outdf <- NULL
  for(i in 1:dim(zp)[1]){
    df <- zp[i, ]
    if(df$z < 2.5){ #non hubs
      if(df$p < 0.05){
        df$roles <- "ultra peripheral"
      }
      else if(df$p < 0.620){
        
        df$roles <- "peripheral"
      }
      else if(df$p < 0.80){
        df$roles <- "non hub connector"
      }
      else{
        df$roles <- "non hub kinless"
      }
    }
    else { # module hubs
      if(df$p < 0.3){
        df$roles <- "provincial hub"
      }
      else if(df$p < 0.75){
        
        df$roles <- "connector hub"
      }
      else {
        df$roles <- "kinless hub"
      }
    }
    if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
  }
  return(outdf)
}

assign_module_roles_4 <- function(zp){
  zp <- na.omit(zp)
  zp$roles <- rep(0, dim(zp)[1])
  outdf <- NULL
  for (i in 1:dim(zp)[1]){
    df = zp[i,]
    if (df$z < 2.5){
      if (df$p < 0.62){
        df$roles = 'Peripherals'
      }
      else{
        df$roles = 'Connectors'
      }
    }
    else{
      if (df$p < 0.62){
        df$roles = 'Module hubs'
      }
      else{
        df$roles = 'Network hubs'
      }
    }
      if(is.null(outdf)){outdf <- df}else{outdf <- rbind(outdf, df)}
    }
    return (outdf)
}
