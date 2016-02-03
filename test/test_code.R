

### A test example run for topographic topic model

nlat <- 10;
nlong <- 20;
topo.config <- expand.grid(1:nlat, 1:nlong); ## create a topological configuration given number of rows and columns

mu_sim <- mu_sim_prior(topo.config,3) ## fix the locations of means/centers for the clusters
lambda_sim <- lambda_sim_prior(K=3) ## fix the spread of the centers for the clusters


initopics_theta <- f_sim_prior(topo.config,3)

omega_iter <- rbind(c(0.5,0.3,0.2),c(0.8,0.1,0.1),c(0.1,0.5,0.4), c(0.2,0.6,0.2), c(0.4,0.3,0.3));

counts <- t(do.call(cbind,lapply(1:dim(omega_iter)[1], function(x) rmultinom(1,1000,prob=omega_iter[x,]%*%t(initopics_theta)))));

library(slam)
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
omega_start <- tpxOmegaStart(X, theta_start)

tmax=10000
wtol=10^(-4)
qn=100
grp=NULL
admix=TRUE
nonzero=FALSE
dcut=-10
verb=1
initopics=NULL
tol=0.1
ord=TRUE

theta <- theta_start

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

## tracking
iter <- 0
dif <- tol+1+qn
update <- TRUE
if(verb>0){
  cat("log posterior increase: " )
  digits <- max(1, -floor(log(tol, base=10))) }

