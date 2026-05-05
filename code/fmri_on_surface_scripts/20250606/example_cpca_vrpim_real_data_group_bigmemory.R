rm(list=ls())
graphics.off()

options(digits.secs=0,
        doFuture.rng.onMisuse='ignore')

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)
library(rsvd)

library(complexlm)
library(pracma)

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



cpca_pim <- function(X_real,X_imag,W_init,niter,beta,bs,verbose=F,prtint=10,W_ref=W_init){

  trmtr <- F
  if(ncol(X_real) > 2*nrow(X_real)){
    trmtr <- T

    dim_ref <- dim(W_init)

    n <- ncol(X_real)
    p <- nrow(X_real)

    W_help <- matrix(NA,nrow=p,ncol=ncol(W_init))
    for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[oi,]),imaginary=as.matrix(X_imag[oi,])),nrow=1,ncol=p) %*% W_init}
    W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
    W_init <- W_help

  }else{
    n <- nrow(X_real)
    p <- ncol(X_real)
  }

  if(length(beta) == 1){beta <- rep(beta,times=niter)}

  R_prev <- diag(1,nrow=ncol(W_init),ncol=ncol(W_init))

  W_prev <- matrix(0,nrow=nrow(W_init),ncol=ncol(W_init))
  W_curr <- W_init

  nb <- round(n/bs)
  b <- (1:n %% nb) + 1
  for(ni in 1:niter){

    b <- b[sample(1:n,size=n)]
    b_list <- list()
    for(ii in 1:nb){
      b_list[[ii]] <- which(b == ii)
    }
    if(any(sapply(b_list,length) > bs)){
      b_rests <- do.call(c,lapply(b_list,function(x) if(length(x) > bs){x[(bs+1):length(x)]}))
      b_rests <- c(b_rests,sample(1:n,size=bs-length(b_rests)))
      b_list <- lapply(b_list,function(x) if(length(x) > bs){x <- x[1:bs];return(x)})
      b_list[[nb+1]] <- b_rests
    }

    A <- matrix(0,nrow=p,ncol=p)
    for(bi in 1:nb){

      if(trmtr){
        X <- matrix(complex(real=as.matrix(t(X_real[,b_list[[bi]]])),imaginary=as.matrix(t(X_imag[,b_list[[bi]]]))),nrow=bs,ncol=p)
      }else{
        X <- matrix(complex(real=as.matrix(X_real[b_list[[bi]],]),imaginary=as.matrix(X_imag[b_list[[bi]],])),nrow=bs,ncol=p)
      }

      X <- X - matrix(rep(apply(X,2,mean),bs),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix

      A <- A + ( Conj(t(X)) %*% ((X)) )/bs
    }

    A <- A/nb

    W_prev_help <- W_curr

    W_curr <- A %*% W_curr - beta[ni]*t(t(inv(R_prev)) %*% t(W_prev))
    R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))
    W_curr <- W_curr %*% inv(R_curr)

    R_prev <- R_curr
    W_prev <- W_prev_help


    if(verbose & ((((ni - 1) %% prtint)) == 0)){
      if(trmtr){
        W_help <- matrix(NA,nrow=dim_ref[1],ncol=dim_ref[2])
        for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[,oi]),imaginary=as.matrix(X_imag[,oi])),nrow=1,ncol=p) %*% W_curr}
        W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
        mccs <- rep(NA,times=ncol(W_ref))
        for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_help[,u]))}
      }else{
        mccs <- rep(NA,times=ncol(W_ref))
        for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
      }
      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Iteration ',ni,' - ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4,width=4),collapse=' - ')))
    }
  }

  if(trmtr){
    W_help <- matrix(NA,nrow=dim_ref[1],ncol=dim_ref[2])
    for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[,oi]),imaginary=as.matrix(X_imag[,oi])),nrow=1,ncol=p) %*% W_curr}
    W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
    return(W_help)
  }else{
    return(W_curr)
  }

}





#
# shuffle_and_adjust_blocks <- function(b,bs,nb){
#   n <- length(b)
#   b <- b[sample(1:n,size=n)]
#   b_list <- list()
#   for(ii in 1:nb){
#     b_list[[ii]] <- which(b == ii)
#   }
#   if(any(sapply(b_list,length) > bs)){
#     b_rests <- do.call(c,lapply(b_list,function(x) if(length(x) > bs){x[(bs+1):length(x)]}))
#     b_rests <- c(b_rests,sample(1:n,size=bs-length(b_rests)))
#     b_list <- lapply(b_list,function(x) if(length(x) > bs){x <- x[1:bs];return(x)})
#     b_list[[nb+1]] <- b_rests
#   }
#   return(b_list)
# }
#
# covariance_per_blocks <- function(X_real,X_imag,trmtr,b_list,bs,nb){
#   if(trmtr){
#     p <- nrow(X_real)
#   }else{
#     p <- ncol(X_real)
#   }
#   A <- matrix(0,nrow=p,ncol=p)
#   for(bi in 1:nb){
#     if(trmtr){
#       X <- matrix(complex(real=as.matrix(t(X_real[,b_list[[bi]]])),imaginary=as.matrix(t(X_imag[,b_list[[bi]]]))),nrow=bs,ncol=p)
#     }else{
#       X <- matrix(complex(real=as.matrix(X_real[b_list[[bi]],]),imaginary=as.matrix(X_imag[b_list[[bi]],])),nrow=bs,ncol=p)
#     }
#     X <- X - matrix(rep(apply(X,2,mean),bs),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
#     A <- A + ( Conj(t(X)) %*% ((X)) )/bs
#   }
#   A <- A/nb
#   return(A)
# }
#
# power_iteration_with_momentum <- function(X_real,X_imag,trmtr,RW_list,beta,b_list,bs,nb){
#
#   R_prev <- RW_list$R
#   W_prev <- RW_list$W_prev
#   W_curr <- RW_list$W_curr
#
#   A <- covariance_per_blocks(X_real=X_real,X_imag=X_imag,trmtr=trmtr,b_list=b_list,bs=bs,nb=nb)
#
#   W_prev_help <- W_curr
#
#   W_curr <- A %*% W_curr - beta*t(t(inv(R_prev)) %*% t(W_prev))
#   R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))
#   W_curr <- W_curr %*% inv(R_curr)
#
#   R_prev <- R_curr
#   W_prev <- W_prev_help
#
#   RW_list <- list(
#                R=R_curr,
#                W_prev=W_prev,
#                W_curr=W_curr
#              )
#
#   return(RW_list)
#
# }
#
# rayleigh_quotient <- function(W,A){
#   mu <- rep(NA,times=ncol(W))
#   for(rqi in 1:length(mu)){
#     mu[rqi] <- (t(W[,rqi]) %*% A %*% W[,rqi])/(t(W[,rqi]) %*% W[,rqi])
#   }
#   mu <- sum(mu)
#   return(mu)
# }
#
# cpca_hb <- function(X_real,X_imag,W_init,niter,bs,verbose=F,W_ref=W_init){
#
#   trmtr <- F
#   if(ncol(X_real) > 2*nrow(X_real)){
#     trmtr <- T
#
#     dim_ref <- dim(W_init)
#
#     n <- ncol(X_real)
#     p <- nrow(X_real)
#
#     W_help <- matrix(NA,nrow=p,ncol=ncol(W_init))
#     for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[oi,]),imaginary=as.matrix(X_imag[oi,])),nrow=1,ncol=p) %*% W_init}
#     W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
#     W_init <- W_help
#
#   }else{
#     n <- nrow(X_real)
#     p <- ncol(X_real)
#   }
#
#   R_prev <- diag(1,nrow=ncol(W_init),ncol=ncol(W_init))
#
#   W_prev <- matrix(0,nrow=nrow(W_init),ncol=ncol(W_init))
#   W_curr <- W_init
#
#   nbf <- round(n/bs)
#   b <- (1:n %% nbf) + 1
#
#   # initial estimate of beta
#   b_list <- shuffle_and_adjust_blocks(b=b,bs=bs,nb=nbf)
#   A <- covariance_per_blocks(X_real=X_real,X_imag=X_imag,trmtr=trmtr,b_list=b_list,bs=bs,nb=nbf)
#
#   mu <- rayleigh_quotient(W_init,A)
#   beta_center <- Mod(mu)^2/4
#   beta_search <- c(2/3,0.99,1.01,1.5)*beta_center
#
#   HB_list <- list()
#   for(hbi in 1:length(beta_search)){
#     HB_list[[hbi]] <- list(RW_list=list(R=R_prev,
#                                         W_prev=W_prev,
#                                         W_curr=W_curr),
#                            beta=beta_search[hbi])
#   }
#
#
#   for(ni in 1:niter){
#
# #     b_list <- list()
# #     for(bi in 1:nb){b_list[[bi]] <- sample(x=1:n,size=bs)}
#
#     b_list <- shuffle_and_adjust_blocks(b=b,bs=bs,nb=nbf)
#
#     for(hbi in 1:length(HB_list)){
#       HB_list[[hbi]]$RW_list <- power_iteration_with_momentum(X_real=X_real,X_imag=X_imag,trmtr=trmtr,RW_list=HB_list[[hbi]]$RW_list,beta=HB_list[[hbi]]$beta,b_list=b_list,bs=bs,nb=nb)
#     }
#
#     if((((ni %% 10)) == 0)){
#
#       mu <- rep(NA,times=length(HB_list))
#       for(hbi in 1:length(HB_list)){mu[hbi] <- rayleigh_quotient(W=HB_list[[hbi]]$RW_list$W_curr,A=A)}
#
#       W_curr <- HB_list[[which.max(mu)]]$RW_list$W_curr
#
#       beta_center <- beta_search[which.max(mu)]
#       beta_search <- c(2/3,0.99,1.01,1.5)*beta_center
#
#       for(hbi in 1:length(HB_list)){HB_list[[hbi]]$beta <- beta_search[hbi]}
#
#     }
#
#     if(verbose & (((ni %% 10)) == 0)){
#       if(trmtr){
#         W_help <- matrix(NA,nrow=dim_ref[1],ncol=dim_ref[2])
#         for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[,oi]),imaginary=as.matrix(X_imag[,oi])),nrow=1,ncol=p) %*% W_curr}
#         W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
#         mccs <- rep(NA,times=ncol(W_ref))
#         for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_help[,u]))}
#       }else{
#         mccs <- rep(NA,times=ncol(W_ref))
#         for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
#       }
#       print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Iteration ',ni,' - ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4,width=4),collapse=' - ')))
#     }
#   }
#
#   if(trmtr){
#     W_help <- matrix(NA,nrow=dim_ref[1],ncol=dim_ref[2])
#     for(oi in 1:n){W_help[oi,] <- matrix(complex(real=as.matrix(X_real[,oi]),imaginary=as.matrix(X_imag[,oi])),nrow=1,ncol=p) %*% W_curr}
#     W_help <- W_help/matrix(rep(apply(W_help,2,function(x) sqrt(sum(x^2))),times=nrow(W_help)),nrow=nrow(W_help),ncol=ncol(W_help),byrow=T)
#     return(W_help)
#   }else{
#     return(W_curr)
#   }
#
# }


cpca_vrpim <- function(X_real,X_imag,W_init,nepoch,niter,beta,bs,nb,bsf,verbose=F,prtint=1,W_ref=W_init,wrkrs=1){

  if(!(length(beta) == 1 | length(beta) == nepoch | length(beta) == nepoch*niter)){stop('length(beta) has to be either 1 or nepoch or nepoch*niter.')}

  if(length(beta) == 1){beta <- rep(beta,times=nepoch*niter)}
  if(length(beta) == nepoch){beta <- do.call(c,lapply(beta,function(x) rep(x,times=niter)))}

  n <- nrow(X_real)
  p <- ncol(X_real)
  k <- ncol(W_init)

  R_prev <- diag(1,nrow=k,ncol=k)

  W_prev <- matrix(0,nrow=p,ncol=k)
  W_curr <- W_init
  W_anch <- W_init

  nbf <- round(n/bsf)
  bf <- (1:n %% nbf) + 1


  for(ne in 1:nepoch){

    bf <- bf[sample(1:n,size=n)]
    bf_list <- list()
    for(ii in 1:nbf){
      bf_list[[ii]] <- which(bf == ii)
    }
    if(any(sapply(bf_list,length) > bsf)){
      bf_rests <- do.call(c,lapply(bf_list,function(x) if(length(x) > bsf){x[(bsf+1):length(x)]}))
      bf_rests <- c(bf_rests,sample(1:n,size=bsf-length(bf_rests)))
      bf_list <- lapply(bf_list,function(x) if(length(x) > bsf){x <- x[1:bsf];return(x)})
      bf_list[[nbf+1]] <- bf_rests
    }

    if(wrkrs == 1){
    A <- matrix(0,nrow=p,ncol=p)
      for(bi in 1:nbf){
        X <- matrix(complex(real=as.matrix(X_real[bf_list[[bi]],]),imaginary=as.matrix(X_imag[bf_list[[bi]],])),nrow=bsf,ncol=p)
        X <- X - matrix(rep(apply(X,2,mean),bsf),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
        A <- A + (Conj(t(X)) %*% ((X)))/bsf
      }
    }else{
      plan(multicore,workers=wrkrs)
      A <- foreach(bi=1:nbf,.combine='+') %dofuture% {
        X <- matrix(complex(real=as.matrix(X_real[bf_list[[bi]],]),imaginary=as.matrix(X_imag[bf_list[[bi]],])),nrow=bsf,ncol=p)
        X <- X - matrix(rep(apply(X,2,mean),bsf),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
        A <- (Conj(t(X)) %*% ((X)))/bsf
        return(A)
      }
      plan(sequential)
    }
    A <- A/nbf

    W_anch <- A %*% W_anch - beta[(ne-1)*niter+1]*t(t(inv(R_prev)) %*% t(W_prev))
    R_anch <- qr.R(qr(rbind(W_anch,W_prev)))
    W_anch <- W_anch %*% inv(R_anch)

    V_anch <- A %*% W_anch

    for(ni in 1:niter){

      b <- sample(x=1:n,size=bs*nb)
      b_list <- list()
      for(bi in 1:nb){
        b_list[[bi]] <- b[((bi-1)*bs+1):(bi*bs)]
      }

      if(wrkrs == 1){
      A <- matrix(0,nrow=p,ncol=p)
        for(bi in 1:nb){
          X <- matrix(complex(real=as.matrix(X_real[b_list[[bi]],]),imaginary=as.matrix(X_imag[b_list[[bi]],])),nrow=bs,ncol=p)
          X <- X - matrix(rep(apply(X,2,mean),bs),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
          A <- A + (Conj(t(X)) %*% ((X)))/bs
        }
      }else{
        plan(multicore,workers=wrkrs)
        A <- foreach(bi=1:nb,.combine='+') %dofuture%{
          X <- matrix(complex(real=as.matrix(X_real[b_list[[bi]],]),imaginary=as.matrix(X_imag[b_list[[bi]],])),nrow=bs,ncol=p)
          X <- X - matrix(rep(apply(X,2,mean),bs),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
          A <- (Conj(t(X)) %*% ((X)))/bs
          return(A)
        }
        plan(sequential)
      }
      A <- A/nb

      W_prev_help <- W_curr

      alpha <- Conj(t(W_curr)) %*% W_anch

      W_curr <- A %*% (W_curr - W_anch %*% alpha) + V_anch %*% alpha - beta[(ne-1)*niter+ni]*Conj(t(Conj(t(inv(R_prev))) %*% Conj(t(W_prev))))
      R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))
      W_curr <- W_curr %*% inv(R_curr)

      R_prev <- R_curr
      W_prev <- W_prev_help

      if(verbose & ((((ni - 1) %% prtint)) == 0)){
        mccs <- rep(NA,times=ncol(W_ref))
        for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Epoch ',formatC(ne,width=1+ceiling(log10(max(c(nepoch,niter))))),' - Iteration ',formatC(ni,width=1+ceiling(log10(max(c(nepoch,niter))))),' - ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4),collapse=' - ')))
      }
    }
  }

  return(W_curr)

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

# parameters for the initial estimate
size_init <- 500 # number of samples for initialization

# # parameters for the pim method
# niter <- 50 # number of iterations
# beta <- 1e-5
# bs <- 20*100 # 50 samples per dataset # block size


path_imag <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/misc/'
plotid <- paste0(sample(c(letters,LETTERS,1:9),size=8),collapse='')
cmap <- rainbow(n_pc)

Sets <- list()
Sets[[1]] <- fold_data[grepl(pattern='^ESO_C',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]
Sets[[2]] <- fold_data[grepl(pattern='^ESO_P',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]
Sets[[3]] <- fold_data

Sizes <- list()

for(set in 1:length(Sets)){

  # parameters for the vrpim method
  nepoch <- 5
  niter <- 7
  beta <- 10^seq(from=-4,to=-10,length.out=nepoch*niter)
  bs <- 100 # block size
  nb <- 10 # number of blocks
  bsf <- length(Sets[[set]])*100 # block size for the anchor estimate

  # (1) prepare the data

  plan(multicore,workers=16)

  for(fold in 1:length(Sets[[set]])){

    if(!dir.exists(paste0(path_work,Sets[[set]][fold]))){dir.create(paste0(path_work,Sets[[set]][fold]),recursive=T)}

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

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Creating the bigmemory data matrices'))

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

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Centering and scaling the data'))

  for(clmn in 1:ncol(data_comb_real)){
    helper <- complex(real=data_comb_real[,clmn],imaginary=data_comb_imag[,clmn])
    mu <- mean(helper)
    sigma2 <- var(helper)
    if(sigma2 > 1e-9){
      helper <- (helper - mu)/sqrt(sigma2)
    }else{
      helper <- (helper - mu)
    }

    data_comb_real[,clmn] <- Re(helper)
    data_comb_imag[,clmn] <- Im(helper)
  }


  blas_set_num_threads(16)
  omp_set_num_threads(16)

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

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating group cPCA using the minibatch power method with momentum'))

  # (4.1) initialization from a small sample of the data

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating initialization on a subset of data'))

  ind_init <- sample(1:data_size[1],size=size_init)
  data_init <- matrix(complex(real=as.matrix(data_comb_real[ind_init,]),imaginary=as.matrix(data_comb_imag[ind_init,])),nrow=size_init,ncol=data_size[2])

  pca_init <- rpca(A=data_init,k=n_pc,center=F,scale=F)
  W <- pca_init$rotation

#   print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating beta value based on the initial pca - beta: ',beta))

  rm('ind_init','data_init')
  gc()

  # keeping the initial values as a reference
  W_init <- W

  # the correlation between the initial estimate and the reference
  mccs_init <- rep(NA,times=n_pc)
  for(u in 1:n_pc){mccs_init[u] <- Mod(cor(x=cpca_opt$rotation[,u],y=W_init[,u]))}

  print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')
  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Modulus of correlation coefficient between initial estimate and reference'))
  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - ',paste0('cPC',1:n_pc,': ',round(mccs_init,digits=4),collapse=' - ')))
  print('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #')

  # (4.2) estimating group cpca by the modified Oja method

#   W <- cpca_pim(X_real=data_comb_real,X_imag=data_comb_imag,W_init=W,niter=niter,beta=beta,bs=bs,verbose=T,prtint=1,W_ref=cpca_opt$rotation)
#   W <- cpca_hb(X_real=data_comb_real,X_imag=data_comb_imag,W_init=W,niter=niter,bs=bs,verbose=T,W_ref=cpca_opt$rotation)
  W <- cpca_vrpim(X_real=data_comb_real,X_imag=data_comb_imag,W_init=W,nepoch=nepoch,niter=niter,beta=beta,bs=bs,nb=nb,bsf=bsf,verbose=F,prtint=1,W_ref=cpca_opt$rotation,wrkrs=16)

  # comparing the results
  mccs_iter  <- rep(NA,times=n_pc)
  for(u in 1:n_pc){mccs_iter[u]  <- Mod(cor(x=cpca_opt$rotation[,u],y=W[,u]))}

  print('### # # #  #  #  #   #   #   #    #    #    #     #     #     #     #    #    #    #   #   #   #  #  #  # # # ###')
  print('Modulus of correlation coefficient between initial estimate and reference')
  print(formatC(mccs_init,format='f',digits=4,width=4))
  print('Modulus of correlation coefficient between current estimate and reference')
  print(formatC(mccs_iter,format='f',digits=4,width=4))
  print('### # # #  #  #  #   #   #   #    #    #    #     #     #     #     #    #    #    #   #   #   #  #  #  # # # ###')

}


