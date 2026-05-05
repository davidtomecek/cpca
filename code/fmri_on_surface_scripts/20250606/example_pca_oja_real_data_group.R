rm(list=ls())
graphics.off()

options(digits.secs=0)

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)
library(rsvd)
#
# get_raut <- function(x,fs=1e3,fsdec=1e3){
#   buff <- 7*fsdec
#   n <- length(x)
#
#   fsdecfilt <- 0.5
#   ord_raut <- buttord(Wp=c(0.005/(fsdecfilt/2), 0.1/(fsdecfilt/2)),Ws=c(0.001/(fsdecfilt/2), 0.14/(fsdecfilt/2)),Rp=3,Rs=20)
#   filt_raut <- butter(n=ord_raut$n,Rs=ord_raut$Rs,w=ord_raut$Wc,type='pass')
#
#   x <- approx(x=seq(from=0,by=1/fs,to=(length(x)-1)/fs),y=x,xout=seq(from=0,by=1/fsdecfilt,to=(length(x)-1)/fsdecfilt),rule=2)$y
#
#   x <- filtfilt(filt=filt_raut$b,
#                 a=filt_raut$a,
#                 x=c(rev(x[2:(buff+1)]),x,rev(x[(length(x)-buff):(length(x)-1)])))[(buff+1):(length(x)+buff)]
#
#   x <- resample(x=x,p=2*fsdec,q=2*fsdecfilt)
#
#   if(length(x) > n){x <- x[1:n]}
#   return(x)
# }
#
# preprocess_dataset <- function(path_data_left,path_data_right,path_mcp,subset_left,subset_right,fs_orig=1,fs_up=3){
#
#   data_left  <- read.nifti1.data(path_data_left)
#   data_right <- read.nifti1.data(path_data_right)
#
#   data_left  <- data_left[subset_left,]
#   data_right <- data_right[subset_right,]
#
#   t_orig <- seq(from=0,by=1/fs_orig,length.out=ncol(data_left))
#   t_up <- seq(from=0,by=1/fs_up,length.out=ncol(data_left)*(fs_up/fs_orig))
#
#   # filtering the data to bandpass 0.005Hz, 0.1Hz
#   crpar <- read.delim(file=path_mcp,header=F,sep='')[,2:7]
#
#   data_left <- foreach(j=1:nrow(data_left),.combine='rbind',.export=c('get_raut')) %dopar% {
# #     m <- stats::lm(data_left[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
#     m <- stats::lm(data_left[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6] + abs(crpar[,1]) + abs(crpar[,2]) + abs(crpar[,3]) + abs(crpar[,4]) + abs(crpar[,5]) + abs(crpar[,6]))
#     x <- resid(m)
#     x <- approx(t_orig,x,t_up)$y
#     x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
#     return(x)
#   }
#
#   data_right <- foreach(j=1:nrow(data_right),.combine='rbind',.export=c('get_raut')) %dopar% {
# #     m <- stats::lm(data_right[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
#     m <- stats::lm(data_right[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6] + abs(crpar[,1]) + abs(crpar[,2]) + abs(crpar[,3]) + abs(crpar[,4]) + abs(crpar[,5]) + abs(crpar[,6]))
#     x <- resid(m)
#     x <- approx(t_orig,x,t_up)$y
#     x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
#     return(x)
#   }
#
#   # adding complex part to the signals using Hilbert transform
#   buff <- 10*fs_up
#   data_left <- foreach(j=1:nrow(data_left),.combine='rbind') %dopar% {
#     x <- c(rev(head(data_left[j,],buff)),data_left[j,],rev(tail(data_left[j,],buff)))
#     x <- hilbert(x=x)
#     x <- x[(buff+1):(ncol(data_left)+buff)]
#     return(x)
#   }
#
#   data_right <- foreach(j=1:nrow(data_right),.combine='rbind') %dopar% {
#     x <- c(rev(head(data_right[j,],buff)),data_right[j,],rev(tail(data_right[j,],buff)))
#     x <- hilbert(x=x)
#     x <- x[(buff+1):(ncol(data_right)+buff)]
#     return(x)
#   }
#
#   return(rbind(data_left,data_right))
#
# }


library(doFuture)

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

  # loading the data
  data_left  <- read.nifti1.data(path_data_left)
  data_right <- read.nifti1.data(path_data_right)

  # focusing on the subset
  data_left  <- data_left[subset_left,]
  data_right <- data_right[subset_right,]

  # auxiliary time variables
  t_orig <- seq(from=0,by=1/fs_orig,length.out=ncol(data_left))
  t_up <- seq(from=0,by=1/fs_up,length.out=ncol(data_left)*(fs_up/fs_orig))

  # motion correction parameters
  crpar <- read.delim(file=path_mcp,header=F,sep='')[,2:7]

  # regressing out the motion correction parameters
  data_left <- foreach(row=1:nrow(data_left),.combine='rbind') %dofuture% {
    m <- stats::lm(data_left[row,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    x <- resid(m)
    return(x)
  }

  data_right <- foreach(row=1:nrow(data_right),.combine='rbind') %dofuture% {
    m <- stats::lm(data_right[row,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    x <- resid(m)
    return(x)
  }

  # normalization of the timesetries by power in the spectrum
  psd_left  <- foreach(row=1:nrow(data_left),.combine=rbind) %dofuture% {spec.pgram(data_left[row,],spans=c(3,5),demean=T,detrend=F,plot=F)$spec}
  psd_right <- foreach(row=1:nrow(data_right),.combine=rbind) %dofuture% {spec.pgram(data_right[row,],spans=c(3,5),demean=T,detrend=F,plot=F)$spec}

  psd_avg <- apply(rbind(psd_left,psd_right),2,mean)
  freq  <- seq(from=0,to=fs_orig/2,length.out=length(psd_avg))

  norm <- sqrt(sum(psd_avg[freq <= 0.15])) # this one looked the best in the visualizations

  rm('psd_left','psd_right','psd_avg')
  gc()

  data_left <- t(apply(data_left,1,function(x) x/norm))
  data_right <- t(apply(data_right,1,function(x) x/norm))

  # upsampling and filtering to the 0.001 - 0.1 Hz band as in raut et al.
  data_left <- foreach(row=1:nrow(data_left),.combine='rbind') %dofuture% {
    x <- approx(x=t_orig,y=data_left[row,],xout=t_up)$y
    x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
    return(x)
  }

  data_right <- foreach(row=1:nrow(data_right),.combine='rbind') %dofuture% {
    x <- approx(x=t_orig,y=data_right[row,],xout=t_up)$y
    x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
    return(x)
  }

  # adding complex part to the signals using Hilbert transform
  buff <- 10*fs_up
  data_left <- foreach(row=1:nrow(data_left),.combine='rbind') %dofuture% {
    x <- c(rev(head(data_left[row,],buff)),data_left[row,],rev(tail(data_left[row,],buff)))
    x <- hilbert(x=x)
    x <- x[(buff+1):(ncol(data_left)+buff)]
    return(x)
  }

  data_right <- foreach(row=1:nrow(data_right),.combine='rbind') %dofuture% {
    x <- c(rev(head(data_right[row,],buff)),data_right[row,],rev(tail(data_right[row,],buff)))
    x <- hilbert(x=x)
    x <- x[(buff+1):(ncol(data_right)+buff)]
    return(x)
  }

  return(rbind(data_left,data_right))

}

library(complexlm)

cpca_oja_cmplx <- function(X,k,R_init=NULL,niter,nu,bs=NULL,ss=NULL,ns=NULL,verbose=F){

  # inputs
  #   X ........ complex data matrix
  #   k ........ number of principal compenents
  #   R_init ... initialization for the rotation matrix
  #   niter .... number of iterations
  #   nu ....... learning rate
  #   bs ....... block size
  #   ss ....... sample size
  #   ns ....... number of samples to draw in the ss case to speed-up the computations
  #
  # outputs
  #   R ... estimated rotation matrix

  if(verbose){prtint <- 10^(floor(log10(niter))-1)}

  trmtr <- F
  if(ncol(X) > 2*nrow(X)){
    trmtr <- T

    if(!is.null(R_init)){
      R_init <- X %*% R_init
      R_init <- R_init/matrix(rep(apply(R_init,2,function(x) sqrt(sum(x^2))),times=nrow(R_init)),nrow=nrow(R_init),ncol=ncol(R_init),byrow=T)
    }

    X <- t(X)
  }

  n <- nrow(X)
  p <- ncol(X)

  if(is.null(R_init)){

    if(verbose){
      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Using random initializaition.'))
    }

    R <- matrix(complex(real=rnorm(n=p*r),imaginary=rnorm(n=p*r)),nrow=p,ncol=r)
    R <- qr.Q(qr(R))

  }else{

    if(nrow(R_init) != ncol(X) | ncol(R_init) != k){
      stop(paste0('The dimension of the initial rotation matrix should be ',ncol(X),'x',k,'.'))
    }

    if(verbose){
      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Using user specified initializaition.'))
    }

    R <- R_init
  }

  if(length(nu) == 1){nu <- rep(nu,times=niter)}
  if(length(nu) != niter){stop('The length of the learning rate vector should be the number of iterations.')}

  if(verbose){
    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Starting estimation'))
    tim <- Sys.time()
  }

  if( (is.null(bs) && is.null(ss)) || (!is.null(bs) && bs == 1) || (!is.null(ss) && ss == 1)){

    for(ni in 1:niter){

      if(is.null(ns)){
        for(no in 1:n){
          A <- X[no,] %*% Conj(t(X[no,]))
          R <- (diag(p) + nu[ni]*A) %*% R
          R <- qr.Q(qr(R))
        }
      }else{

        if(ns > n){stop('The size of subsample is larger than number of observations.')}


        si <- sample(x=1:n,size=ns)
        for(no in 1:ns){
          A <- X[si[no],] %*% Conj(t(X[si[no],]))
          R <- (diag(p) + nu[ni]*A) %*% R
          R <- qr.Q(qr(R))
        }


      }

      if(verbose && (ni-1) %% prtint == 0){
        difftim <- Sys.time() - tim
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
      }
      if(verbose && round(ni == prtint/4 | ni == prtint/2)){
        difftim <- Sys.time() - tim
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
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
          print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }
        if(verbose && round(ni == prtint/4 | ni == prtint/2)){
          difftim <- Sys.time() - tim
          print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }

      }
    }

    if(!is.null(ss)){
      for(ni in 1:niter){
        si <- sample(1:n,size=ss)
        A <- (Conj(t(X[si,])) %*% ((X[si,])))/ss
        R <- (diag(p) + nu[ni]*A) %*% R
        R <- qr.Q(qr(R))

        if(verbose && (ni-1) %% prtint == 0){
          difftim <- Sys.time() - tim
          print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }
        if(verbose && round(ni == prtint/4 | ni == prtint/2)){
          difftim <- Sys.time() - tim
          print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Iteration ',ni,', elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
        }

      }
    }
  }

  if(verbose){
    difftim <- Sys.time() - tim
    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Finished, elapsed time ',as.numeric(round(difftim,digits=2)),' ',attr(difftim,'units')))
  }

  if(trmtr){
    R <- X %*% R
    R <- R/matrix(rep(apply(R,2,function(x) sqrt(sum(x^2))),times=nrow(R)),nrow=nrow(R),ncol=ncol(R),byrow=T)
  }

  return(R)

}

path_data <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/rest/'
fold_data <- list.dirs(path_data,full.names=F,recursive=F)
fold_data <- fold_data[grep(pattern='^ESO_',x=fold_data)]

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
path_save <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

dec <- 'fsaverage5'
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

# force preprocessing even though there are some preprocessed data?
force_preproc <- F

# parameters of data preprocessing
fs_orig <- 1
fs_up <- 2
fs <- fs_up

# number of estimated principal components
n_pc <- 3

# parameters for the Oja method
ninit <- 1.6e3 # number of samples for initialization
niter <- 4 # number of iterations
nepoch <- 128 # number of epochs

# # logarithmically decreasing learning rate
# nu_start <- 1e-2 # learning rate - starting value
# nu_end <- 1e-5 # learning rate - final value
# nu_logdiff <- (log10(nu_start) - log10(nu_end))/nepoch # learning rate - step on log scale
# nu <- nu_start

# # learning rate that is small, gets bigger and then decreases, to make use of the initialization
# nu_low <- 1e-2
# nu_high <- 5e-1
# nu <- (-cos(seq(from=0,to=2*pi,length.out=nepoch))/2 + 0.5)*(nu_high - nu_low) + nu_low

# a different learning rate

nu_nodes <- matrix(
              c( -5, 16.0,
                  1, 8.0,
                 10, 4.0,
                 90, 2.0,
                120, 1.0,
                140, 0.5),
            ncol=2,byrow=T)
nu_nodes <- as.data.frame(nu_nodes)
names(nu_nodes) <- c('x','y')
splinterp <- splinefun(x=nu_nodes$x,y=log(nu_nodes$y),method='fmm')
nu <- exp(splinterp(1:nepoch))

# bs <- 4e2 # block size
ss <- 256 # sample size
# ns <- 32

ndatinit <- 3 # number of datasets for initialization
ndat <- 15 # number of datasets per epoch
nbatch_max <- 32 # number of randomly selected combinations

pi <- 8 # printing something at everty pi-th iteration

path_imag <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/misc/'
plotid <- paste0(sample(c(letters,LETTERS,1:9),size=8),collapse='')
cmap <- rainbow(n_pc)

Sets <- list()
Sets[[1]] <- fold_data[grepl(pattern='^ESO_C',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]

for(set in 1:length(Sets)){


  # (1) prepare the data

#   library(doParallel)
#   compclust <- makeForkCluster(4)
#   registerDoParallel(compclust)

  plan(multicore,workers=4)

  for(fold in 1:length(Sets[[set]])){

    if(!dir.exists(paste0(path_work,Sets[[set]][fold]))){dir.create(dir.exists(paste0(path_work,Sets[[set]][fold])),recursive=T)}

    if(file.exists(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData')) & !force_preproc){

      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Using already preprocessed data for ',Sets[[set]][fold]))

    }else{

      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Preprocessing data for ',Sets[[set]][fold]))

      # preprocessing the data
      data_preproc <- preprocess_dataset(path_data_left=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.sm0.',dec,'.lh.nii.gz'),
                                         path_data_right=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.sm0.',dec,'.rh.nii.gz'),
                                         path_mcp=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.mcdat'),
                                         subset_left=veri_left,
                                         subset_right=veri_right,
                                         fs_orig=fs_orig,
                                         fs_up=fs_up)
      # saving the preprocessed data
      save(data_preproc,file=paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))

      # deleting the preprocessed data
      rm('data_preproc')
      gc()

    }
  }

#   registerDoSEQ()
#   stopCluster(compclust)

  plan(sequential)

  blas_set_num_threads(12)
  omp_set_num_threads(12)

  # (2) estimate the reference solution using rpca from the rsvd package

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating reference group cPCA using rpca from package rsvd'))

  data_comb <- foreach(fold=1:length(Sets[[set]]),.combine='rbind') %do% {
    load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
    data_preproc <- t(data_preproc)
    return(data_preproc)
  }

  cpca_opt <- rpca(A=data_comb,k=n_pc,center=F,scale=F)

  rm('data_comb')
  gc()


  # (3) simulate the Oja method in group setting

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating group cPCA using the modified Oja method'))

  ####
  # in order to speed up the experiments, I'm not going to do the loading data during every batch, in this smaller data set it should be that memory taxing
  Data_preproc <- foreach(fold=1:length(Sets[[set]])) %do% {
    load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
    return(data_preproc)
  }
  ####

  # initialization from a small sample of the data

  data_comb <- foreach(fold=sample(x=1:length(Sets[[set]]),size=ndatinit),.combine='rbind') %do% {
    ####
    #     load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
    ####
    data_preproc <- Data_preproc[[fold]]
    data_preproc <- t(data_preproc)
    return(data_preproc)
  }

  pca_init <- rpca(A=data_comb[sample(1:nrow(data_comb),size=ninit),],k=n_pc,center=F,scale=F)
  R <- pca_init$rotation

  # keeping the initial values as a reference
  R_init <- R

  # preparing the indices for the dataset combinations
  datind <- combn(x=1:length(Sets[[set]]),m=ndat)

  # logging the correlations between the rsvd and moja solutions per each iteration
  Mccs_moja <- list()
  Mccs_moja_help <- list()
  moja <- 1

  mccs_init <- rep(NA,times=n_pc)
  for(u in 1:n_pc){mccs_init[u] <- Mod(cor(x=cpca_opt$rotation[,u],y=R_init[,u]))}

  mccs_moja_bsf <- rep(0,times=n_pc)

  # estimating group cpca by the modified Oja method
  for(epoch in 1:nepoch){

#     print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Starting epoch ',epoch,' out of ',nepoch))

    # reshuffling the data sets so that the order is different in each epoch
    batch_mix <- sample(x=1:ncol(datind),size=min(c(nbatch_max,ncol(datind))))

    for(batch in 1:min(c(nbatch_max,ncol(datind)))){

      if(batch == 1 | batch == min(c(nbatch_max,ncol(datind))) | (batch - 1) %% pi == 0){
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Epoch ',formatC(epoch,width=ceiling(log10(nepoch))),'/',nepoch,' - batch ',formatC(batch,width=ceiling(log10(min(c(nbatch_max,ncol(datind)))))),'/',min(c(nbatch_max,ncol(datind))),' - learning rate ',formatC(nu[epoch],format='g')))
      }

      data_comb <- foreach(fold=datind[,batch_mix[batch]],.combine='rbind') %do% {
        ####
        #     load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
        ####
        data_preproc <- Data_preproc[[fold]]
        data_preproc <- t(data_preproc)
        return(data_preproc)
      }

      R <- cpca_oja_cmplx(X=data_comb,k=n_pc,R_init=R,niter=niter,nu=nu[epoch],ss=ss,verbose=F)

      ####
      Mccs_moja_help[[moja]] <- Mod(foreach(u=1:n_pc,.combine='c') %do% {cor(x=cpca_opt$rotation[,u],y=R[,u])})
      mccs_moja_bsf <- apply(rbind(mccs_moja_bsf,Mccs_moja_help[[moja]]),2,max)
      print(round(Mccs_moja_help[[moja]],digits=3))
      moja <- moja + 1
      ####

    }

#     # updating the learning rate that changes with epochs
#     nu <- 10^(log10(nu) - nu_logdiff)

    # comparing the results
    mccs_moja  <- rep(NA,times=n_pc)
    for(u in 1:n_pc){mccs_moja[u]  <- Mod(cor(x=cpca_opt$rotation[,u],y=R[,u]))}
    Mccs_moja[[epoch]] <- mccs_moja

    print('### # # #  #  #  #   #   #   #    #    #    #     #     #     #     #    #    #    #   #   #   #  #  #  # # # ###')
    print('Modulus of correlation coefficient between initial estimate and reference')
    print(round(mccs_init,digits=3))
    print('Modulus of correlation coefficient between current moja estimate and reference')
    print(round(mccs_moja,digits=3))
    print('Best-so-far correlation coefficient')
    print(round(mccs_moja_bsf,digits=3))
    print('### # # #  #  #  #   #   #   #    #    #    #     #     #     #     #    #    #    #   #   #   #  #  #  # # # ###')

    png(filename=paste0(path_imag,'mccs_moja_',plotid,'.png'),width=5+(2*epoch*min(c(nbatch_max,ncol(datind)))/1000),height=5,res=300,units='in')
    plot(x=1:(moja - 1),xlab='Epoch number',
         y=sapply(Mccs_moja_help,function(x) x[1]),ylim=c(0,1),ylab='MCC mOja x rSVD',
         col=cmap[1],lwd=0.6,
         type='l',frame=F,xaxt='n')
    axis(side=1,at=(1:epoch)*min(c(nbatch_max,ncol(datind))),labels=1:epoch)
    abline(h=(0:10)/10,lty=3)
    for(pc in 2:n_pc){
      lines(x=1:(moja - 1),
            y=sapply(Mccs_moja_help,function(x) x[pc]),
            col=cmap[pc],lwd=0.6)
    }
    legend('topright',legend=paste0('cPC',1:n_pc),lwd=1.4,col=cmap,bty='n')
    dev.off()

  }

  Mccs_moja <- do.call(rbind,Mccs_moja)

  cmap <- rainbow(n_pc)
  png(filename=paste0(path_imag,'mccs_moja_',gsub(pattern='-',replacement='',x=Sys.Date()),'.png'),width=9,height=5,res=300,units='in')
  plot(x=1:nepoch,xlab='Epoch number',
       y=Mccs_moja[,1],ylim=c(0,1),ylab='MCC mOja x rSVD',
       col=cmap[1],lwd=1.4,
       type='l',frame=F)
  for(pc in 2:n_pc){
    lines(x=1:nepoch,
          y=Mccs_moja[,pc],
          col=cmap[pc],lwd=1.4)
  }
  legend('topright',legend=paste0('cPC',1:n_pc),lwd=1.4,col=cmap,bty='n')
  dev.off()
}
