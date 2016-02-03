##### Estimation for Topic Models ######
## intended main function; provides defaults and fits topic model for the user defined K
topo.topics <- function(counts, K, shape=NULL, initopics=NULL, z_options=1, tol=0.1,
                  ord=TRUE, verb=1, reflect=TRUE, ...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL, admix=TRUE, nonzero=FALSE, dcut=-10
{
  ## Convert the counts matrix to triplet matrix
  X <- CheckCounts(counts)
  p <- ncol(X)
  n <- nrow(X)

  if(verb > 0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the list of candidate K values
  if(K <=1){ stop(cat("use K values >= 2")) }
  cat(sprintf("\nFitting a topographical topic model with %d topics \n", K))

  ## Null model log probability
  sx <- sum(X)
  qnull <- col_sums(X)/sx
  null <- sum( X$v*log(qnull[X$j]) ) - 0.5*(n+p)*(log(sx) - log(2*pi))


  ## initialize
  f_start <- f_sim_prior(topo.config,K);
  theta_start <- t(do.call(rbind, lapply(1:K, function(k) gtools::rdirichlet(1, f_start[,k]))));
  
  ## initialize

  fit <- topo.tpxfit(counts=counts, X=X, theta=theta_start, f=f_start, z_options, tol=tol, verb=verb, 
                admix=admix, grp=grp, tmax=tmax, wtol=wtol, qn=qn);


  #initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1], shape, verb)

  ## either search for marginal MAP K and return bayes factors, or just fit
  ## tpx <- tpxSelect(X, K, bf, initopics, alpha=shape, tol, kill, verb, ...)
  ## K <- tpx$K

  ## clean up and out
  if(ord){ worder <- order(col_sums(fit$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(fit$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(fit$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  f = matrix(fit$f[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  mu = matrix(fit$mu[worder,], nrow=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  lambda = matrix(fit$lambda[worder,], nrow=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }
  
  ## topic object
  out <- list(K=K, theta=theta, omega=omega, f=f, mu=mu, lambda=lambda, loglik=fit$L, X=X, null=null)
  class(out) <- "topics"
  invisible(out) }


