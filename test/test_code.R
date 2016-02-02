

### A test example run for topographic topic model

nlat <- 10;
nlong <- 20;
topo.config <- expand.grid(1:nlat, 1:nlong); ## create a topological configuration given number of rows and columns

mu_sim <- mu_sim_prior(topo.config,3) ## fix the locations of means/centers for the clusters
lambda_sim <- lambda_sim_prior(K=3) ## fix the spread of the centers for the clusters


initopics_theta <- f_sim_prior(topo.config,3)

omega_iter <- rbind(c(0.5,0.3,0.2),c(0.8,0.1,0.1),c(0.1,0.5,0.4));

counts <- t(do.call(cbind,lapply(1:dim(omega_iter)[1], function(x) rmultinom(1,1000,prob=omega_iter[x,]%*%t(initopics_theta)))));
