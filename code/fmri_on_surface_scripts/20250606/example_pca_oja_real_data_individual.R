rm(list=ls())
graphics.off()

options(digits.secs=0)

library(RhpcBLASctl)
blas_set_num_threads(4)
omp_set_num_threads(4)

library(freesurferformats)
library(gsignal)
library(rsvd)

library(doParallel)
compclust <- makeForkCluster(4)
registerDoParallel(compclust)

# library(doMC)
# registerDoMC(8)

get_raut <- function(x,fs=1e3,fsdec=1e3){
  buff <- 7*fsdec
  n <- length(x)

  fsdecfilt <- 0.5
  ord_raut <- buttord(Wp=c(0.005/(fsdecfilt/2), 0.1/(fsdecfilt/2)),Ws=c(0.001/(fsdecfilt/2), 0.14/(fsdecfilt/2)),Rp=3,Rs=20)
  filt_raut <- butter(n=ord_raut$n,Rs=ord_raut$Rs,w=ord_raut$Wc,type='pass')

  x <- approx(x=seq(from=0,by=1/fs,to=(length(x)-1)/fs),y=x,xout=seq(from=0,by=1/fsdecfilt,to=(length(x)-1)/fsdecfilt),rule=2)$y

  x <- filtfilt(filt=filt_raut$b,
                a=filt_raut$a,
                x=c(rev(x[2:(buff+1)]),x,rev(x[(length(x)-buff):(length(x)-1)])))[(buff+1):(length(x)+buff)]

  x <- resample(x=x,p=2*fsdec,q=2*fsdecfilt)

  if(length(x) > n){x <- x[1:n]}
  return(x)
}

preprocess_dataset <- function(path_data_left,path_data_right,path_mcp,subset_left,subset_right,fs_orig=1,fs_up=3){

  data_left  <- read.nifti1.data(path_data_left)
  data_right <- read.nifti1.data(path_data_right)

  data_left  <- data_left[subset_left,]
  data_right <- data_right[subset_right,]

  t_orig <- seq(from=0,by=1/fs_orig,length.out=ncol(data_left))
  t_up <- seq(from=0,by=1/fs_up,length.out=ncol(data_left)*(fs_up/fs_orig))

  # filtering the data to bandpass 0.005Hz, 0.1Hz
  crpar <- read.delim(file=path_mcp,header=F,sep='')[,2:7]

  data_left <- foreach(j=1:nrow(data_left),.combine='rbind',.export=c('get_raut')) %dopar% {
#     m <- stats::lm(data_left[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    m <- stats::lm(data_left[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6] + abs(crpar[,1]) + abs(crpar[,2]) + abs(crpar[,3]) + abs(crpar[,4]) + abs(crpar[,5]) + abs(crpar[,6]))
    x <- resid(m)
    x <- approx(t_orig,x,t_up)$y
    x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
    return(x)
  }

  data_right <- foreach(j=1:nrow(data_right),.combine='rbind',.export=c('get_raut')) %dopar% {
#     m <- stats::lm(data_right[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    m <- stats::lm(data_right[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6] + abs(crpar[,1]) + abs(crpar[,2]) + abs(crpar[,3]) + abs(crpar[,4]) + abs(crpar[,5]) + abs(crpar[,6]))
    x <- resid(m)
    x <- approx(t_orig,x,t_up)$y
    x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
    return(x)
  }

  # adding complex part to the signals using Hilbert transform
  buff <- 10*fs_up
  data_left <- foreach(j=1:nrow(data_left),.combine='rbind') %dopar% {
    x <- c(rev(head(data_left[j,],buff)),data_left[j,],rev(tail(data_left[j,],buff)))
    x <- hilbert(x=x)
    x <- x[(buff+1):(ncol(data_left)+buff)]
    return(x)
  }

  data_right <- foreach(j=1:nrow(data_right),.combine='rbind') %dopar% {
    x <- c(rev(head(data_right[j,],buff)),data_right[j,],rev(tail(data_right[j,],buff)))
    x <- hilbert(x=x)
    x <- x[(buff+1):(ncol(data_right)+buff)]
    return(x)
  }

  return(rbind(data_left,data_right))

}


library(complexlm)

cpca_oja_cmplx <- function(X,k,R_init=NULL,niter,nu,bs=NULL,ss=NULL,verbose=F){

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

  if(verbose){prtint <- 10^(floor(log10(niter))-1)}

  n <- nrow(X)
  p <- ncol(X)

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

  if(verbose){
    print(paste0(Sys.time(),': Starting estimation'))
    tim <- Sys.time()
  }

  if( (is.null(bs) && is.null(ss)) || (!is.null(bs) && bs == 1) || (!is.null(ss) && ss == 1)){

    for(ni in 1:niter){
      for(no in 1:n){
        A <- X[no,] %*% Conj(t(X[no,]))
        R <- (diag(p) + nu[ni]*A) %*% R
        R <- qr.Q(qr(R))
      }

      if(verbose && (ni-1) %% prtint == 0){
        difftim <- Sys.time() - tim
        print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
      }
      if(verbose && round(ni == prtint/4 | ni == prtint/2)){
        difftim <- Sys.time() - tim
        print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
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

        if(verbose && (ni-1) %% prtint == 0){
          difftim <- Sys.time() - tim
          print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }
        if(verbose && round(ni == prtint/4 | ni == prtint/2)){
          difftim <- Sys.time() - tim
          print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }

      }
    }

    if(!is.null(ss)){
      for(ni in 1:niter){
        si <- sample(1:n,size=ss)
        A <- (Conj(t(X[si])) %*% ((X[si,])))/ss
        R <- (diag(p) + nu[ni]*A) %*% R
        R <- qr.Q(qr(R))

        if(verbose && (ni-1) %% prtint == 0){
          difftim <- Sys.time() - tim
          print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }
        if(verbose && round(ni == prtint/4 | ni == prtint/2)){
          difftim <- Sys.time() - tim
          print(paste0(Sys.time(),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }

      }
    }
  }

  if(verbose){
    difftim <- Sys.time() - tim
    print(paste0(Sys.time(),': Finished, elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
  }

  return(R)

}

path_data <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/rest/'
fold_data <- list.dirs(path_data,full.names=F,recursive=F)
fold_data <- fold_data[grep(pattern='^ESO_',x=fold_data)]

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
path_save <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

dec <- 'fsaverage6'
path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')
path_labe <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/label/')

atlas_left  <- read.fs.annot(paste0(path_labe,'lh.aparc.annot'))
atlas_right <- read.fs.annot(paste0(path_labe,'rh.aparc.annot'))

hemi_left  <- 1:length(atlas_left$vertices)
hemi_right <- (1:length(atlas_left$vertices)) + length(atlas_left$vertices)

veri_left  <-  which(!(atlas_left$label_names %in% c('unknown','corpuscallosum')))
veri_right <- which(!(atlas_right$label_names %in% c('unknown','corpuscallosum')))

rm('atlas_left','atlas_right')
gc()

# parameters of data preprocessing
fs_orig <- 1
fs_up <- 2
fs <- fs_up

# number of estimated principal components
n_pc <- 5

# parameters for the Oja method
ninit <- 2048 # number of samples for initialization
niter <- 128 # number of iterations
nu <- 10^seq(from=-2,to=-5,length.out=niter) # learning rate
bs <- 2048 # block size
nsim <- 32 # number of simulations

for(i in sample(x=1:length(fold_data),size=1)){

  print(paste0(Sys.time(),': Preprocessing data for ',fold_data[i]))

#   # loading fMRI data
#   data_left <- read.nifti1.data(paste0(path_data,fold_data[i],'/rest/001/fmcpr.sm0.',dec,'.lh.nii.gz'))
#   data_right <- read.nifti1.data(paste0(path_data,fold_data[i],'/rest/001/fmcpr.sm0.',dec,'.rh.nii.gz'))
#
#   # loading motion correction parameters
#   crpar <- read.delim(file=paste0(path_data,fold_data[i],'/rest/001/fmcpr.mcdat'),header=F,sep='')[,2:7]

  # preprocessing the data
  data_preproc <- preprocess_dataset(path_data_left=paste0(path_data,fold_data[i],'/rest/001/fmcpr.sm0.',dec,'.lh.nii.gz'),
                                     path_data_right=paste0(path_data,fold_data[i],'/rest/001/fmcpr.sm0.',dec,'.rh.nii.gz'),
                                     path_mcp=paste0(path_data,fold_data[i],'/rest/001/fmcpr.mcdat'),
                                     subset_left=veri_left,
                                     subset_right=veri_right,
                                     fs_orig=fs_orig,
                                     fs_up=fs_up)

  # complex pca - reference solution by rsvd
  print(paste0(Sys.time(),': Estimating cPCA for ',fold_data[i]))
  cpca_opt <- rpca(A=data_preproc,k=n_pc,center=F,scale=F)

  # complex pca - experiments with the Oja method
  print(paste0(Sys.time(),': Oja cPCA experiment for ',fold_data[i]))

  registerDoSEQ()
  stopCluster(compclust)

  Res <- foreach(ns=1:nsim) %do% {

    print(paste0(Sys.time(),': Starting simulation ',ns))

    # initialization from a small sample of the data
    pca_init <- rpca(A=data_preproc[sample(1:nrow(data_preproc),size=ninit),],k=n_pc,center=F,scale=F)
    R <- pca_init$rotation

    # keeping the initial values as a reference
    R_init <- R

    # estimating cpca by Oja method
    R <- cpca_oja_cmplx(X=data_preproc,k=n_pc,R_init=R_init,niter=niter,nu=nu,bs=bs,verbose=T)

    # comparing the results
    mccs_init <- rep(NA,times=n_pc)
    mccs_oja  <- rep(NA,times=n_pc)

    for(u in 1:n_pc){
      mccs_init[u] <- Mod(cor(x=cpca_opt$rotation[,u],y=R_init[,u]))
      mccs_oja[u]  <- Mod(cor(x=cpca_opt$rotation[,u],y=R[,u]))
    }

    return(list(mccs_init=mccs_init,mccs_oja=mccs_oja))
  }

  mccs_init <- do.call(rbind,lapply(Res,function(x) x$mccs_init))
  mccs_oja  <- do.call(rbind,lapply(Res,function(x) x$mccs_oja))


  # number of initial solutions with correlation coefficient greater the 0.95
  print(
    apply(mccs_init,2,function(x) sum(x > 0.95))
  )
  # proportion of initial solutions with correlation coefficient greater the 0.95
  print(
    apply(mccs_init,2,function(x) sum(x > 0.95))/nsim
  )

  # number of oja solutions with correlation coefficient greater the 0.95
  print(
    apply(mccs_oja,2,function(x) sum(x > 0.95))
  )
  # proportion of oja solutions with correlation coefficient greater the 0.95
  print(
    apply(mccs_oja,2,function(x) sum(x > 0.95))/nsim
  )

  # number of initial solutions with correlation coefficient greater the 0.975
  print(
    apply(mccs_init,2,function(x) sum(x > 0.975))
  )
  # proportion of initial solutions with correlation coefficient greater the 0.975
  print(
    apply(mccs_init,2,function(x) sum(x > 0.975))/nsim
  )

  # number of oja solutions with correlation coefficient greater the 0.975
  print(
    apply(mccs_oja,2,function(x) sum(x > 0.975))
  )
  # proportion of oja solutions with correlation coefficient greater the 0.975
  print(
    apply(mccs_oja,2,function(x) sum(x > 0.975))/nsim
  )

  # # saving the results
  # save(list=c('fs','hemi_left','hemi_right','veri_left','veri_right','cpca_opt'),file=paste0(path_save,fold_data[i],'/',dec,'_cpca_opt.RData'))

}

