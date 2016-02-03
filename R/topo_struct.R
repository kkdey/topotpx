
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

CheckCounts <- function(fcounts){
  if(class(fcounts)[1] == "TermDocumentMatrix"){ fcounts <- t(fcounts) }
  if(is.null(dimnames(fcounts)[[1]])){ dimnames(fcounts)[[1]] <- paste("doc",1:nrow(fcounts)) }
  if(is.null(dimnames(fcounts)[[2]])){ dimnames(fcounts)[[2]] <- paste("wrd",1:ncol(fcounts)) }
  empty <- row_sums(fcounts) == 0
  if(sum(empty) != 0){
    fcounts <- fcounts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(fcounts))
}


topo.tpxOmegaStart <- function(X, theta)
{
  if(!inherits(X,"simple_triplet_matrix")){ stop("X needs to be a simple_triplet_matrix.") }
  omega <- try(tcrossprod_simple_triplet_matrix(X, solve(t(theta)%*%theta)%*%t(theta)), silent=TRUE )
  if(inherits(omega,"try-error")){ return( matrix( 1/ncol(theta), nrow=nrow(X), ncol=ncol(theta) ) ) }
  omega[omega <= 0] <- .5
  return( normalize(omega, byrow=TRUE) )
}

topo.tpxlpost <- function(counts, omega, theta, f)
{
  prob <- omega%*%t(theta);
  K <- ncol(theta);
  omega[omega<=0] <- 1e-20;
  theta[theta<=0] <- 1e-20;
  prob[prob <=0] <- 1e-20;
  L <- sum(log(omega))/K + sum(f*log(theta))+sum(counts*log(prob));
  return(L)
}

z_construct <- function(counts, omega_iter, theta_iter, z_options=c(1,2))
{
  if(z_options==2){
    row_total <- rowSums(counts);
    z_est <- round(sweep(theta_iter, MARGIN=1, colSums(sweep(omega_iter, MARGIN = 1, row_total, "*")), "*"));
  }
  if(z_options==1){
    z_est <- (t(omega_iter) %*% (counts/(omega_iter %*% theta_iter)))*theta_iter;
  }
  return(z_est)
}

## ** called from topics.R (predict) and tpx.R
## Conditional solution for topic weights given theta
topo.tpxweights <- function(n, p, xvo, wrd, doc, start, theta, verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000)
{
  K <- ncol(theta)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start) 
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              X = as.double(xvo),
              theta = as.double(theta),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="maptpx")
  return(t(matrix(omega$W, nrow=ncol(theta), ncol=n))) }


topo.tpxfit <- function(counts, X, theta, f, z_options, tol, verb,
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
  omega <- topo.tpxOmegaStart(X,theta);
  
  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb>0){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }
  
  Y <- NULL # only used for qn > 0 
  Q0 <- col_sums(X)/sum(X)
  L <- topo.tpxlpost(counts, omega=omega, theta=theta, f=f) 
  # if(is.infinite(L)){ L <- sum( (log(Q0)*col_sums(X))[Q0>0] ) }
  
  ## Iterate towards MAP
  while( update  && iter < tmax ){
    
    ## sequential quadratic programming for conditional Y solution
    
    if(admix && wtol > 0){ Wfit <- topo.tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc,
                                              start=omega, theta=theta,  verb=0, nef=TRUE, wtol=wtol, tmax=50) }
    else{ Wfit <- omega }
    
    z_construct <- z_construct(counts, omega_iter = Wfit, theta_iter = t(theta), z_options = 1);
    
    theta_fit <- t(normalize(z_construct + t(f)));
    
    move <- list(theta=theta_fit, omega=Wfit);
    
    QNup <- topo.tpxQN(move=move, counts=counts, Y=Y, X=X, f=f, verb=verb, admix=admix, grp=grp, doqn=qn-dif)
    move <- QNup$move
    Y <- QNup$Y
    
    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
      
      z_construct <- z_construct(counts, omega_iter = omega, theta_iter = t(theta), z_options = 1);
      theta_fit <- t(normalize(z_construct + t(f)));
      move <- list(theta=theta_fit, omega=Wfit);
      QNup$L <-  topo.tpxlpost(counts, move$omega, move$theta, f);
    }
    
    ## calculate dif
    dif <- (QNup$L-L)
    
    L <- QNup$L
    
    
    ## check convergence
    if(abs(dif) < tol){
      if(sum(abs(theta-move$theta)) < tol){ update = FALSE } }
    
    ## print
    if(verb>0 && (iter-1)%%ceiling(10/verb)==0 && iter>0){
      cat( paste( round(dif,digits), #" (", sum(abs(theta-move$theta)),")",
                  ", ", sep="") ) }
    
    ## heartbeat for long jobs
    if(((iter+1)%%1000)==0){ 
      cat(sprintf("p %d iter %d diff %g\n",
                  nrow(theta), iter+1,round(dif))) }
    
    iter <- iter+1
    theta <- move$theta;
    omega <- move$omega;
    
    mu_fit <- do.call(rbind,lapply(1:K, function(k) colSums(sweep(topo.config, MARGIN=1, theta[,k], "*"))));
    
    lambda_fit <- do.call(rbind, lapply(1:K, function(k) colSums(sweep((topo.config-sweep(cbind(rep(1,dim(topo.config)[1]),rep(1,dim(topo.config)[1])),MARGIN = 2, (mu_fit[k,]), "*"))^2, MARGIN = 1, theta[,k],"*"))))
    
    f <- do.call(cbind,lapply(1:K, function(k)
    {
      loglik_temp <- sapply(1:dim(topo.config)[1], 
                            function(s) dmvnorm(as.numeric(topo.config[s,]),
                                                as.numeric(mu_fit[k,]), diag(lambda_fit[k], dim(topo.config)[2]), 
                                                log=TRUE));
      weights <- exp(loglik_temp - max(loglik_temp));
      prop_weights <- weights/ sum(weights);
      return(prop_weights)
    }))
    
  }
  
  L <- topo.tpxlpost(counts, omega=omega, theta=theta, f=f)
  
  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }
  
  out <- list(theta=theta, omega=omega, K=K, f=f, L=L, mu=mu_fit, lambda=lambda_fit, iter=iter)
  invisible(out) }


topo.tpxQN <- function(move, counts, Y, X, f, verb, admix, grp, doqn)
{
  ## always check likelihood
  L <- topo.tpxlpost(counts, omega=move$omega, theta=move$theta, f=f); 
  
  if(doqn < 0){ return(list(move=move, L=L, Y=Y)) }
  
  ## update Y accounting
  Y <- cbind(Y, topo.tpxToNEF(theta=move$theta, omega=move$omega))
  if(ncol(Y) < 3){ return(list(Y=Y, move=move, L=L)) }
  if(ncol(Y) > 3){ warning("mis-specification in quasi-newton update; please report this bug.") }
  
  ## Check quasinewton secant conditions and solve F(x) - x = 0.
  U <- as.matrix(Y[,2]-Y[,1])
  V <- as.matrix(Y[,3]-Y[,2])
  sUU <- sum(U^2)
  sVU <- sum(V*U)
  Ynew <- Y[,3] + V*(sVU/(sUU-sVU)) 
  qnup <- topo.tpxFromNEF(Ynew, n=nrow(move$omega),
                     p=nrow(move$theta), K=ncol(move$theta))
  
  ## check for a likelihood improvement
  Lqnup <- try(topo.tpxlpost(counts, omega=qnup$omega, theta=qnup$theta, f=f), silent=TRUE)
  
  if(inherits(Lqnup, "try-error")){
    if(verb>10){ cat("(QN: try error) ") }
    return(list(Y=Y[,-1], move=move, L=L)) }
  
  if(verb>10){ cat(paste("(QN diff ", round(Lqnup-L,3), ")\n", sep="")) }
  
  if(Lqnup < L){
    return(list(Y=Y[,-1], move=move, L=L)) }
  else{
    L <- Lqnup
    Y <- cbind(Y[,2],Ynew)
    return( list(Y=Y, move=qnup, L=L) )
  }
}

## fast computation of sparse P(X) for X>0
topo.tpxQ <- function(theta, omega, doc, wrd){
  
  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(theta)){stop("theta/omega mis-match in tpxQ") }
  
  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(theta)),
            K = as.integer(ncol(theta)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            theta = as.double(theta),
            q = double(length(wrd)),
            PACKAGE="topotpx" )
  
  return( out$q ) }

## model and component likelihoods for mixture model
topo.tpxMixQ <- function(X, omega, theta, grp=NULL, qhat=FALSE){
  if(is.null(grp)){ grp <- rep(1, nrow(X)) }
  K <- ncol(omega)
  n <- nrow(X)
  mixhat <- .C("RmixQ",
               n = as.integer(nrow(X)),
               p = as.integer(ncol(X)),
               K = as.integer(K),
               N = as.integer(length(X$v)),
               B = as.integer(nrow(omega)),
               cnt = as.double(X$v),
               doc = as.integer(X$i-1),
               wrd = as.integer(X$j-1),
               grp = as.integer(as.numeric(grp)-1),
               omega = as.double(omega),
               theta = as.double(theta),
               Q = double(K*n),
               PACKAGE="topotpx")
  ## model and component likelihoods
  lQ <- matrix(mixhat$Q, ncol=K)
  lqlhd <- log(row_sums(exp(lQ)))
  lqlhd[is.infinite(lqlhd)] <- -600 # remove infs
  if(qhat){
    qhat <- exp(lQ-lqlhd)
    ## deal with numerical overload
    infq <- row_sums(qhat) < .999
    if(sum(infq)>0){
      qhat[infq,] <- 0
      qhat[n*(apply(matrix(lQ[infq,],ncol=K),1,which.max)-1) + (1:n)[infq]] <- 1 }
  }
  return(list(lQ=lQ, lqlhd=lqlhd, qhat=qhat)) }


## functions to move theta/omega to and from NEF.  
topo.tpxToNEF <- function(theta, omega){
  n <- nrow(omega)
  p <- nrow(theta)
  K <- ncol(omega)
  return(.C("RtoNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=double((p-1)*K + n*(K-1)),
            theta=as.double(theta), tomega=as.double(t(omega)),
            PACKAGE="topotpx")$Y)
}

## 'From' NEF representation back to probabilities
topo.tpxFromNEF <- function(Y, n, p, K){
  bck <- .C("RfromNEF",
            n=as.integer(n), p=as.integer(p), K=as.integer(K),
            Y=as.double(Y), theta=double(K*p), tomega=double(K*n),
            PACKAGE="topotpx")
  return(list(omega=t( matrix(bck$tomega, nrow=K) ), theta=matrix(bck$theta, ncol=K)))
}
