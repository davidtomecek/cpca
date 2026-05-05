rm(list=ls())
graphics.off()

library(foreach)

library(doParallel)
clst <- makeForkCluster(4)
registerDoParallel(clst)

library(complexlm)

cpca_oja_cmplx <- function(X,k,R_init=NULL,niter,nu,bs=NULL,ss=NULL){
  
  # inputs
  #   X ........ complex data matrix
  #   k ........ number of principal compenents
  #   R_init ... initialization for the rotation matrix
  #   niter .... number of iterations
  #   nu ....... learning rate
  #   bs ....... block size
  #   ss ....... sample size
  #
  # outputs
  #   R ... estimated rotation matrix
  
  n <- nrow(X)
  
  if(is.null(R)){
    R <- matrix(complex(real=rnorm(n=p*r),imaginary=rnorm(n=p*r)),nrow=p,ncol=r)
    R <- qr.Q(qr(R))
  }else{
    if(nrow(R_init) != ncol(X) | ncol(R_init) != k){
      stop(paste0('The dimension of the initial rotation matrix should be ',ncol(X),'x',k,'.'))
    }
  }
  
  if(length(nu) == 1){nu <- rep(nu,times=niter)}
  if(length(nu) != niter){stop('The length of the learning rate vector should be the number of iterations.')}
  
  if( (is.null(bs) && is.null(ss)) || (!is.null(bs) && bs == 1) || (!is.null(ss) && ss == 1)){
    
    for(ni in 1:niter){
      for(no in 1:n){
        A <- X[no,] %*% Conj(t(X[no,]))
        R <- (diag(p) + nu[ni]*A) %*% R
        R <- qr.Q(qr(R))
      }
    }
    
  }else{
    
    if(!is.null(bs)){
      nb <- round(n/bs)
      b <- (1:n %% nb) + 1
      for(ni in 1:niter){
        b <- b[sample(1:n,size=n)]
        for(bi in 1:nb){
          A <- (Conj(t(X[b == bi,])) %*% ((X[b == bi,])))/sum(b == bi)
          R <- (diag(p) + nu[ni]*A) %*% R
          R <- qr.Q(qr(R))
        }
      }
    }
    
    if(!is.null(ss)){
      for(ni in 1:niter){
        si <- sample(1:n,size=ss)
        A <- (Conj(t(X[si])) %*% ((X[si,])))/ss
        R <- (diag(p) + nu[ni]*A) %*% R
        R <- qr.Q(qr(R))
      }
    }
  }
  
  return(R)
  
}

# assessing the Oja method for complex PCA convergence in a simulation

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'

n <- 1024 # number of observations
p <- 128 # number of variables
r <- 5 # number of components
s <- 256 # number of samples for initialization
niter <- 2e3 # number of iterations
nu <- 10^seq(from=-2,to=-5,length.out=niter) # learning rate
bs <- 32 # block size
nsim <- 100 # number of simulations

# mccs_init <- matrix(NA,nrow=nsim,ncol=r)
# mccs_oja <- matrix(NA,nrow=nsim,ncol=r)

# for(ns in 1:nsim){

Res <- foreach(ns=1:nsim) %dopar% {

  print(paste0(Sys.time(),': Starting iteration ',ns))
  
  # random data matrix
  D <- matrix(complex(real=rnorm(n=n*p),imaginary=rnorm(n=n*p)),nrow=n,ncol=p)
  D <- scale(D,center=T,scale=T)
  
  # reference pca solution
  pca_ref <- prcomp(D,center=F,scale=F)

  # # initializing the first loading vector and normalizing it to unit length
  # R <- matrix(complex(real=runif(n=p*r,min=0,max=1),imaginary=runif(n=p*r,min=0,max=1)),nrow=p,ncol=r)
  # R <- qr.Q(qr(R))
  
  # trying to get a better initialization
  pca_init <- svd(D[sample(1:n,size=s),])
  R <- pca_init$v[,1:r]

  # keeping the initial values for plotting
  R_init <- R

  # estimating cpca by Oja method
  R <- cpca_oja_cmplx(X=D,k=r,R_init=R_init,niter=niter,nu=nu,bs=bs)

  # # comparing the results
  # for(u in 1:r){
  #   mccs_init[ns,u] <- Mod(cor(x=pca_ref$rotation[,u],y=R_init[,u]))
  #   mccs_oja[ns,u] <- Mod(cor(x=pca_ref$rotation[,u],y=R[,u]))
  # }

  mccs_init <- rep(NA,times=r)
  mccs_oja  <- rep(NA,times=r)
  
  # comparing the results
  for(u in 1:r){
    mccs_init[u] <- Mod(cor(x=pca_ref$rotation[,u],y=R_init[,u]))
    mccs_oja[u]  <- Mod(cor(x=pca_ref$rotation[,u],y=R[,u]))
  }
  
  print(paste0(Sys.time(),': Finished iteration ',ns))
  
  return(list(mccs_init=mccs_init,mccs_oja=mccs_oja))
}

mccs_init <- do.call(rbind,lapply(Res,function(x) x$mccs_init))
mccs_oja  <- do.call(rbind,lapply(Res,function(x) x$mccs_oja))

# some summaries

# number of initial solutions with correlation coefficient greater the 0.95
apply(mccs_init,2,function(x) sum(x > 0.95))
# proportion of initial solutions with correlation coefficient greater the 0.95
apply(mccs_init,2,function(x) sum(x > 0.95))/nsim

# number of oja solutions with correlation coefficient greater the 0.95
apply(mccs_oja,2,function(x) sum(x > 0.95))
# proportion of oja solutions with correlation coefficient greater the 0.95
apply(mccs_oja,2,function(x) sum(x > 0.95))/nsim

save.image(paste0(path_work,'oja_simulation_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))

# what's next
# find a measure of the goodness of fit for pca
# guess how many runs of the algorithm are necessary to get to the correct solution
# try working with the actual data from resting state fmri experiments 

stopCluster(clst)