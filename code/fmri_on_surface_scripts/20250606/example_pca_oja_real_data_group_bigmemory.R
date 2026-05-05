rm(list=ls())
graphics.off()

options(digits.secs=0)

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)
library(rsvd)

library(bigmemory)


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

cpca_expe <- function(X,k,R_init,niter,nu,beta,bs){

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

  if(length(nu) == 1){nu <- rep(nu,times=niter)}
  if(length(beta) == 1){beta <- rep(beta,times=niter)}

  R_prev <- R_init
  R_curr <- R_init

  nb <- round(n/bs)
  b <- (1:n %% nb) + 1
  for(ni in 1:niter){
    b <- b[sample(1:n,size=n)]

    A <- matrix(0,nrow=p,ncol=p)
    for(bi in 1:nb){
      A <- A + (Conj(t(X[b == bi,])) %*% ((X[b == bi,])))/sum(b == bi)
    }

    A <- A/nb

    R_prev_help <- R_curr

    R_curr <- nu[ni]*A %*% R_curr - beta[ni]*R_prev
    R <- qr.Q(qr(R))

    R_prev <- R_prev_help

  }

  if(trmtr){
    R_curr <- X %*% R_curr
    R_curr <- R_curr/matrix(rep(apply(R_curr,2,function(x) sqrt(sum(x^2))),times=nrow(R_curr)),nrow=nrow(R_curr),ncol=ncol(R_curr),byrow=T)
  }

  return(R_curr)

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
size_init <- 500 # number of samples for initialization
niter <- 1 # number of iterations
nepoch <- 64 # number of epochs
nbatch <- 30 # number of batches
# size_batch <- 200 # size of batch for training, use this when using some of the simple sampling strategies

# parameters for sampling data in blocks
size_block <- 50 # size of block of data per subset
overlap <- T # whether the blocks of data can overlap

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
#
# nu_nodes <- matrix(
#               c(  -16, 32.000,
#                    16, 16.000,
#                    32,  8.000,
#                    64,  4.000,
#                   128,  2.000,
#                   256,  1.000,
#                   512,  0.500,
#                  1024,  0.250,
#                  2048,  0.125),
#             ncol=2,byrow=T)
# nu_nodes <- matrix(
#               c(  -5, 1e+0,
#                    1, 1e+0,
#                  100, 5e-1,
#                  150, 1e-1,
#                  200, 5e-2,
#                  300, 1e-2),
#             ncol=2,byrow=T)
# nu_nodes <- matrix(
#               c(  -5, 1e-3,
#                    1, 1e-3,
#                  100, 5e-4,
#                  150, 1e-4,
#                  200, 5e-5,
#                  300, 1e-5),
#             ncol=2,byrow=T)
# nu_nodes <- matrix(
#               c(  -5, 1e-18,
#                    1, 1e-18,
#                  100, 5e-19,
#                  150, 1e-19,
#                  200, 5e-20,
#                  300, 1e-20),
#             ncol=2,byrow=T)
# nu_nodes <- as.data.frame(nu_nodes)
# names(nu_nodes) <- c('x','y')
# splinterp <- splinefun(x=nu_nodes$x,y=log(nu_nodes$y),method='fmm')
# nu <- exp(splinterp(1:nepoch))

# learning rate & momentum in cpca_expe
nu <- 10^seq(from=-3,to=-5,length.out=nepoch)
beta <- 10^seq(from=-1,to=-3,length.out=nepoch)

bs <- 1000 # block size
# ss <- 2048 # sample size
# ns <- 32

pi <- 5 # printing something at everty pi-th iteration

path_imag <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/misc/'
plotid <- paste0(sample(c(letters,LETTERS,1:9),size=8),collapse='')
cmap <- rainbow(n_pc)

Sets <- list()
Sets[[1]] <- fold_data[grepl(pattern='^ESO_C',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]

Sizes <- list()

for(set in 1:length(Sets)){


  # (1) prepare the data

  plan(multicore,workers=4)

  for(fold in 1:length(Sets[[set]])){

    if(!dir.exists(paste0(path_work,Sets[[set]][fold]))){dir.create(dir.exists(paste0(path_work,Sets[[set]][fold])),recursive=T)}

    if(file.exists(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData')) & !force_preproc){

      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Using already preprocessed data for ',Sets[[set]][fold]))

      load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
      Sizes[[fold]] <- dim(data_preproc)
      rm('data_preproc')
      gc()

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

      # logging data sizes for the construction of the bigmemory matrix
      Sizes[[fold]] <- dim(data_preproc)

      # saving the preprocessed data
      save(data_preproc,file=paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))

      # deleting the preprocessed data
      rm('data_preproc')
      gc()

    }
  }


  plan(sequential)

  # (2) construct the bigmemory matrices for real and imaginary parts

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Creating the bigmemory data matrices.'))

  data_size <- c(sum(sapply(Sizes,function(x) x[2])),Sizes[[1]][1])

  if(length(list.files(path_work,pattern='^data_comb_')) > 0){unlink(paste0(path_work,list.files(path_work,pattern='^data_comb_')))}

  data_comb_real <- filebacked.big.matrix(nrow=data_size[1],ncol=data_size[2],backingpath=path_work,backingfile='data_comb_real',descriptorfile='data_comb_real.desc')
  data_comb_imag <- filebacked.big.matrix(nrow=data_size[1],ncol=data_size[2],backingpath=path_work,backingfile='data_comb_imag',descriptorfile='data_comb_imag.desc')

  curr_row <- 1
  for(fold in 1:length(Sets[[set]])){

    load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))

    dimnames(data_preproc) <- NULL

    data_comb_real[curr_row:(curr_row+Sizes[[set]][2]-1),] <- Re(t(data_preproc))
    data_comb_imag[curr_row:(curr_row+Sizes[[set]][2]-1),] <- Im(t(data_preproc))

    curr_row <- curr_row + Sizes[[set]][2]
  }

  blas_set_num_threads(12)
  omp_set_num_threads(12)

  # (3) estimate the reference solution using rpca from the rsvd package

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating reference group cPCA using rpca from package rsvd'))
#
#   data_comb <- foreach(fold=1:length(Sets[[set]]),.combine='rbind') %do% {
#     load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
#     data_preproc <- t(data_preproc)
#     return(data_preproc)
#   }

  data_comb <- matrix(complex(real=as.matrix(data_comb_real),imaginary=as.matrix(data_comb_imag)),nrow=data_size[1],ncol=data_size[2])

  cpca_opt <- rpca(A=data_comb,k=n_pc,center=F,scale=F)

  rm('data_comb')
  gc()


  # (4) simulate the Oja method in group setting

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating group cPCA using the modified Oja method'))

  # (4.1) initialization from a small sample of the data

  ind_init <- sample(1:data_size[1],size=size_init)
  data_init <- matrix(complex(real=as.matrix(data_comb_real[ind_init,]),imaginary=as.matrix(data_comb_imag[ind_init,])),nrow=size_init,ncol=data_size[2])

  pca_init <- rpca(A=data_init,k=n_pc,center=F,scale=F)
  R <- pca_init$rotation

  rm('ind_init','data_init')
  gc()

  # keeping the initial values as a reference
  R_init <- R

  # logging the correlations between the rsvd and moja solutions per each iteration
  Mccs_moja <- list()
  Mccs_moja_help <- list()
  moja <- 1

  mccs_init <- rep(NA,times=n_pc)
  for(u in 1:n_pc){mccs_init[u] <- Mod(cor(x=cpca_opt$rotation[,u],y=R_init[,u]))}

  mccs_moja_bsf <- rep(0,times=n_pc)

  # indices for mixed sampling of blocks for each subset
  ind_batch_cs <- cumsum(sapply(Sizes,function(x) x[2]))
  ind_batch_cs <- c(0,ind_batch_cs[1:(length(ind_batch_cs)-1)])

  ind_batch_helper <- list()
  if(overlap){
    for(ind in 1:length(Sizes)){
      ind_block <- floor(seq(from=1,to=Sizes[[ind]][2]-size_block+1,length.out=nbatch))
      helper <- list()
      for(bind in 1:nbatch){
        helper[[bind]] <- ind_block[bind]:(ind_block[bind]+size_block-1) + ind_batch_cs[ind]
      }
      ind_batch_helper[[ind]] <- helper
    }
  }else{
    for(ind in 1:length(Sizes)){
      helper <- split(x=1:Sizes[[ind]][2] + ind_batch_cs[ind],f=((1:Sizes[[ind]][2])-1) %/% size_block)
      ind_batch_helper[[ind]] <- helper
    }
  }

  # (4.2) estimating group cpca by the modified Oja method
  for(epoch in 1:nepoch){

#     print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Starting epoch ',epoch,' out of ',nepoch))

    # to ensure that all the data are presented to the learning algorithm during the optimization
    ind_batch_helper <- lapply(ind_batch_helper,function(x) sample(x))

    for(batch in 1:nbatch){

      if(batch == 1 | batch == nbatch | (batch - 1) %% pi == 0){
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Epoch ',formatC(epoch,width=ceiling(log10(nepoch))),'/',nepoch,' - batch ',formatC(batch,width=ceiling(log10(nbatch))),'/',nbatch,' - learning rate ',formatC(nu[epoch],format='g')))
      }

#       data_comb <- foreach(fold=datind[,batch],.combine='rbind') %do% {
#         ####
#         #     load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
#         ####
#         data_preproc <- Data_preproc[[fold]]
#         data_preproc <- t(data_preproc)
#         return(data_preproc)
#       }

#       # selecting randomly from all data
#       ind_batch <- sample(1:data_size[1],size=size_batch)

#       # selecting certain number of samples from each dataset
#       ind_batch_size <- rep(size_batch %/% length(Sizes),times=length(Sizes))
#       ind_batch_rem <- size_batch %% length(Sizes)
#       for(ind in sample(1:length(Sizes),size=ind_batch_rem)){ind_batch_size[ind] <- ind_batch_size[ind] + 1}
#
#       ind_batch <- list()
#       ind_batch_cs <- cumsum(sapply(Sizes,function(x) x[2]))
#       ind_batch_cs <- c(0,ind_batch_cs[1:(length(ind_batch_cs)-1)])
#       for(ind in 1:length(Sizes)){ind_batch[[ind]] <- sort(sample(1:Sizes[[ind]][2],size=ind_batch_size[ind]) + ind_batch_cs[ind])}
#
#       ind_batch <- do.call(c,ind_batch)

#       # selecting blocks of data in order
#       ind_batch_size <- size_batch %/% length(Sizes)
#       ind_batch <- list()
#       ind_batch_cs <- cumsum(sapply(Sizes,function(x) x[2]))
#       ind_batch_cs <- c(0,ind_batch_cs[1:(length(ind_batch_cs)-1)])
#       for(ind in 1:length(Sizes)){ind_batch[[ind]] <- (1:ind_batch_size) + (batch-1)*ind_batch_size + ind_batch_cs[ind]}
#       ind_batch <- do.call(c,ind_batch)

#       # selecting the samples and creating the complex data matrix from the real and imaginary parts
#       data_batch <- matrix(complex(real=as.matrix(data_comb_real[ind_batch,]),imaginary=as.matrix(data_comb_imag[ind_batch,])),nrow=size_batch,ncol=data_size[2])

      # selecting continuous blocks of data from each subset at random
      ind_batch <- do.call(c,lapply(ind_batch_helper,function(x) x[[batch]]))

      # selecting the samples and creating the complex data matrix from the real and imaginary parts
      data_batch <- matrix(complex(real=as.matrix(data_comb_real[ind_batch,]),imaginary=as.matrix(data_comb_imag[ind_batch,])),nrow=length(ind_batch),ncol=data_size[2])


#       R <- cpca_oja_cmplx(X=data_batch,k=n_pc,R_init=R,niter=niter,nu=nu[epoch],bs=bs,verbose=F)
      R <- cpca_expe(X=data_batch,k=n_pc,R_init=R,niter=niter,nu=nu[epoch],beta=beta[epoch],bs=bs)

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

    png(filename=paste0(path_imag,'mccs_moja_',plotid,'.png'),width=5+(2*epoch*nbatch/1000),height=5,res=300,units='in')
    plot(x=1:(moja - 1),xlab='Epoch number',
         y=sapply(Mccs_moja_help,function(x) x[1]),ylim=c(0,1),ylab='MCC mOja x rSVD',
         col=cmap[1],lwd=0.6,
         type='l',frame=F,xaxt='n')
    axis(side=1,at=(1:epoch)*nbatch,labels=1:epoch)
    abline(h=(0:10)/10,lty=3)
    for(pc in 2:n_pc){
      lines(x=1:(moja - 1),
            y=sapply(Mccs_moja_help,function(x) x[pc]),
            col=cmap[pc],lwd=0.6)
    }
    legend('topright',legend=paste0('cPC',1:n_pc),lwd=1.4,col=cmap,bty='n')
    dev.off()

  }

  rm('ind_batch','data_batch')
  gc()

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
