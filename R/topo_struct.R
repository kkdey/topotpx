
### Topographical Topic models

## We assume all the cells are equally likely to be the mean centers for the topogaphical maps

mu_sim_prior <- function(topo.config, K)
{
  ## topo.config is a V \times D matrix, with V cells and for each cell, specifying D dimensions (latitude, longitude, elevation etc)
  return(topo.config[sample(1:dim(topo.config)[1], K, replace=FALSE),])
}

lambda_sim_prior <- function(scale=0.5,K)
{
  ## scale represents the parameter for the Inverse-gamma distribution from which we simulate
  return(1/(rgamma(K,scale,scale)))
}

## simulating the spatial proportional intensities 

f_sim_prior <- function(topo.config, K, scale=0.5)
{
  mu_sim <- mu_sim_prior(topo.config,K);
  lambda_sim <- lambda_sim_prior(scale,K);
  initopics_theta <- do.call(cbind,lapply(1:K, function(k)
              {
                    loglik_temp <- sapply(1:dim(topo.config)[1], 
                                        function(s) dmvnorm(as.numeric(topo.config[s,]),
                                          as.numeric(mu_sim[k,]), diag(lambda_sim[k], dim(topo.config)[2]), 
                                              log=TRUE));
                    weights <- exp(loglik_temp - max(loglik_temp));
                    prop_weights <- weights/ sum(weights);
                    return(prop_weights)
  }))
  
  return(initopics_theta)
}


tpxOmegaStart <- function(X, theta)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
  if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
  omega[omega <= 0] <- .5
  return( normalize(omega, byrow=TRUE) )
}

topotpx::tpxlpost <- function(counts, omega, theta, f)
{
  K <- ncol(theta)
  L <- 0
  L <- L+sum(log(omega))/K + sum(f*log(theta));
  return(L)
}

tootpx::tpxfit <- function(X, theta, alpha, tol, verb,
                   admix, grp, tmax, wtol, qn)
{
  ## inputs and dimensions
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix") }
  K <- ncol(theta)
  n <- nrow(X)
  p <- ncol(X)
  m <- row_sums(X)
  
  ## recycle these in tpcweights to save time
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  
  ## Initialize
  omega <- tpxweights(n=n, p=p, xvo=xvo, wrd=wrd, doc=doc, start=tpxOmegaStart(X,theta), theta=theta)
  if(!admix){ omega <- matrix(apply(omega,2, function(w) tapply(w,grp,mean)), ncol=K) }
  
  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb>0){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }
  
  Y <- NULL # only used for qn > 0 
  Q0 <- col_sums(X)/sum(X)
  L <- tpxlpost(counts, omega=omega, theta=theta, f=f) 
  # if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }
  
  ## Iterate towards MAP
  while( update  && iter < tmax ){
    
    ## sequential quadratic programming for conditional Y solution
    
    if(admix && wtol > 0){ Wfit <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                              start=omega, theta=theta,  verb=0, nef=TRUE, wtol=wtol, tmax=20) }
    else{ Wfit <- omega }
    
    