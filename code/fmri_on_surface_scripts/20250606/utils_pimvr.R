library(gsignal)

library(complexlm)
library(pracma)

library(doFuture)
library(future)
library(listenv)



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
















cpca_pim <- function(X_real,X_imag,W_init,niter,beta,bs,nb,onr=1e2,verbose=F,prtint=1,W_ref=W_init,wrkrs=1){

  if(!(length(beta) == 1 | length(beta) == niter)){stop('length(beta) has to be either 1 or niter.')}

  if(length(beta) == 1){beta <- rep(beta,times=niter)}

  n <- nrow(X_real)
  p <- ncol(X_real)
  k <- ncol(W_init)

  R_prev <- diag(1,nrow=k,ncol=k)

  W_prev <- matrix(0,nrow=p,ncol=k)
  W_curr <- W_init

  for(ni in 1:niter){

    b <- sample(x=1:n,size=bs*nb)
    b_list <- list()
    for(bi in 1:nb){
      b_list[[bi]] <- b[((bi-1)*bs+1):(bi*bs)]
    }

    W_prev_help <- W_curr

    hlpr_w <- listenv()
    rb <- c(seq(from=1,to=p,by=onr),p+1)
    for(ri in 1:(length(rb)-1)){
      hlpr_w[[ri]] %<-% {
        A_sub <- matrix(0,nrow=rb[ri+1]-rb[ri]+1,ncol=p)
        hlpr_b <- listenv()
        for(bi in 1:nb){
          hlpr_b[[bi]] %<-% {
            x1 <- matrix(
                    complex(real=     as.vector(X_real[b_list[[bi]],rb[ri]:(rb[ri+1]-1)]),
                            imaginary=as.vector(X_imag[b_list[[bi]],rb[ri]:(rb[ri+1]-1)])),
                    nrow=bs,ncol=rb[ri+1]-rb[ri])
            x2 <- matrix(
                    complex(real=     as.vector(X_real[b_list[[bi]],]),
                            imaginary=as.vector(X_imag[b_list[[bi]],])),
                    nrow=bs,ncol=p)
            Conj(t(x1)) %*% x2
          }
        }
        A_sub <- Reduce('+',as.list(hlpr_b))
        A_sub <- A_sub/nb
        A_sub %*% W_curr
      }
    }
    W_curr <- do.call(rbind,as.list(hlpr_w))

    W_curr <- W_curr - beta[ni]*Conj(t(Conj(t(inv(R_prev))) %*% Conj(t(W_prev))))
    R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))
    W_curr <- W_curr %*% inv(R_curr)

    R_prev <- R_curr
    W_prev <- W_prev_help

    if(verbose & ((((ni - 1) %% prtint)) == 0)){
      if(!is.null(W_ref)){
        mccs <- rep(NA,times=ncol(W_ref))
        for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Iteration ',formatC(ni,width=1+ceiling(log10(niter))),' - ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4),collapse=' - ')))
      }else{
        print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Iteration ',formatC(ni,width=1+ceiling(log10(niter)))))
      }
    }
  }

  plan(sequential)

  return(W_curr)

}
















# version where the whole covariance matrix is estimated

# cpca_pimvr <- function(X_real,X_imag,W_init,nepoch,niter,beta,bs,nb,bsf,vbs=1e2,verbose=F,prtint=1,W_ref=W_init,wrkrs=1){
#
#   if(!(length(beta) == 1 | length(beta) == nepoch | length(beta) == nepoch*niter)){stop('length(beta) has to be either 1 or nepoch or nepoch*niter.')}
#
#   if(length(beta) == 1){beta <- rep(beta,times=nepoch*niter)}
#   if(length(beta) == nepoch){beta <- do.call(c,lapply(beta,function(x) rep(x,times=niter)))}
#
#   n <- nrow(X_real)
#   p <- ncol(X_real)
#   k <- ncol(W_init)
#
#   R_prev <- diag(1,nrow=k,ncol=k)
#
#   W_prev <- matrix(0,nrow=p,ncol=k)
#   W_curr <- W_init
#   W_anch <- W_init
#
#   nbf <- round(n/bsf)
#   bf <- (1:n %% nbf) + 1
#
#
#   for(ne in 1:nepoch){
#
#     bf <- bf[sample(1:n,size=n)]
#     bf_list <- list()
#     for(ii in 1:nbf){
#       bf_list[[ii]] <- which(bf == ii)
#     }
#     if(any(sapply(bf_list,length) > bsf)){
#       bf_rests <- do.call(c,lapply(bf_list,function(x) if(length(x) > bsf){x[(bsf+1):length(x)]}))
#       bf_rests <- c(bf_rests,sample(1:n,size=bsf-length(bf_rests)))
#       bf_list <- lapply(bf_list,function(x) if(length(x) > bsf){x <- x[1:bsf];return(x)})
#       bf_list[[nbf+1]] <- bf_rests
#     }
#
#     if(wrkrs == 1){
#       A <- matrix(0,nrow=p,ncol=p)
#       for(bi in 1:nbf){
#         X <- matrix(complex(real=as.matrix(X_real[bf_list[[bi]],]),imaginary=as.matrix(X_imag[bf_list[[bi]],])),nrow=bsf,ncol=p)
#         X <- X - matrix(rep(apply(X,2,mean),bsf),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
#         A <- A + (Conj(t(X)) %*% ((X)))/bsf
#       }
#     }else{
#       plan(list(sequential,tweak(multicore,workers=wrkrs)))
#       A <- matrix(0,nrow=p,ncol=p)
#       for(bi in 1:nbf){
#         hlpr1 <- listenv()
#         for(vi1 in seq(from=1,to=p,by=vbs)){
#           hlpr1[[vi1]] %<-% {
#             hlpr2 <- listenv()
#             X1 <- matrix(
#                     complex(real     =as.matrix(X_real[bf_list[[bi]],vi1:min(c(vi1+vbs-1,p))]),
#                             imaginary=as.matrix(X_imag[bf_list[[bi]],vi1:min(c(vi1+vbs-1,p))])),
#                     nrow=bsf,ncol=min(c(vbs,p-vi1+1)))
#             X1 <- apply(X1,2,function(x) x - mean(x))
#             for(vi2 in seq(from=1,to=p,by=vbs)){
#               hlpr2[[vi2]] %<-% {
#                 if(vi2 < vi1){
#                   matrix(0,nrow=min(c(vbs,p-vi1+1)),ncol=min(c(vbs,p-vi2+1)))
#                 }else{
#                   X2 <- matrix(
#                           complex(real     =as.matrix(X_real[bf_list[[bi]],vi2:min(c(vi2+vbs-1,p))]),
#                                   imaginary=as.matrix(X_imag[bf_list[[bi]],vi2:min(c(vi2+vbs-1,p))])),
#                           nrow=bsf,ncol=min(c(vbs,p-vi2+1)))
#                   X2 <- apply(X2,2,function(x) x - mean(x))
#                   (Conj(t(X1)) %*% X2)/bsf
#                 }
#               }
#             }
#             do.call(cbind,as.list(hlpr2))
#           }
#         }
#         A <- A + do.call(rbind,as.list(hlpr1))
#       }
#       plan(sequential)
#     }
#     A <- A/nbf
#     A[lower.tri(A)] <- Conj(t(A))[lower.tri(A)]
#
#     W_anch <- A %*% W_anch - beta[(ne-1)*niter+1]*t(t(inv(R_prev)) %*% t(W_prev))
#     R_anch <- qr.R(qr(rbind(W_anch,W_prev)))
#     W_anch <- W_anch %*% inv(R_anch)
#
#     V_anch <- A %*% W_anch
#
#     for(ni in 1:niter){
#
#       b <- sample(x=1:n,size=bs*nb)
#       b_list <- list()
#       for(bi in 1:nb){
#         b_list[[bi]] <- b[((bi-1)*bs+1):(bi*bs)]
#       }
#
#       if(wrkrs == 1){
#         A <- matrix(0,nrow=p,ncol=p)
#         for(bi in 1:nb){
#           X <- matrix(complex(real=as.matrix(X_real[b_list[[bi]],]),imaginary=as.matrix(X_imag[b_list[[bi]],])),nrow=bs,ncol=p)
#           X <- X - matrix(rep(apply(X,2,mean),bs),ncol=p,byrow=T) # centering the data, necessary for the estimation of the covariance matrix
#           A <- A + (Conj(t(X)) %*% ((X)))/bs
#         }
#       }else{
#
#         plan(list(sequential,tweak(multicore,workers=wrkrs)))
#         A <- matrix(0,nrow=p,ncol=p)
#         for(bi in 1:nb){
#           hlpr1 <- listenv()
#           for(vi1 in seq(from=1,to=p,by=vbs)){
#             hlpr1[[vi1]] %<-% {
#               hlpr2 <- listenv()
#               X1 <- matrix(
#                       complex(real     =as.matrix(X_real[b_list[[bi]],vi1:min(c(vi1+vbs-1,p))]),
#                               imaginary=as.matrix(X_imag[b_list[[bi]],vi1:min(c(vi1+vbs-1,p))])),
#                       nrow=bs,ncol=min(c(vbs,p-vi1+1)))
#               X1 <- apply(X1,2,function(x) x - mean(x))
#               for(vi2 in seq(from=1,to=p,by=vbs)){
#                 hlpr2[[vi2]] %<-% {
#                   if(vi2 < vi1){
#                     matrix(0,nrow=min(c(vbs,p-vi1+1)),ncol=min(c(vbs,p-vi2+1)))
#                   }else{
#                     X2 <- matrix(
#                             complex(real     =as.matrix(X_real[b_list[[bi]],vi2:min(c(vi2 + vbs-1,p))]),
#                                     imaginary=as.matrix(X_imag[b_list[[bi]],vi2:min(c(vi2 + vbs-1,p))])),
#                             nrow=bs,ncol=min(c(vbs,p-vi2+1)))
#                     X2 <- apply(X2,2,function(x) x - mean(x))
#                     (Conj(t(X1)) %*% X2)/bs
#                   }
#                 }
#               }
#               do.call(cbind,as.list(hlpr2))
#             }
#           }
#           A <- A + do.call(rbind,as.list(hlpr1))
#         }
#         plan(sequential)
#       }
#       A <- A/nb
#       A[lower.tri(A)] <- Conj(t(A))[lower.tri(A)]
#
#       W_prev_help <- W_curr
#
#       alpha <- Conj(t(W_curr)) %*% W_anch
#
#       W_curr <- A %*% (W_curr - W_anch %*% alpha) + V_anch %*% alpha - beta[(ne-1)*niter+ni]*Conj(t(Conj(t(inv(R_prev))) %*% Conj(t(W_prev))))
#       R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))
#       W_curr <- W_curr %*% inv(R_curr)
#
#       R_prev <- R_curr
#       W_prev <- W_prev_help
#
#       if(verbose & ((((ni - 1) %% prtint)) == 0)){
#         if(!is.null(W_ref)){
#           mccs <- rep(NA,times=ncol(W_ref))
#           for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
#           print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Epoch ',formatC(ne,width=1+ceiling(log10(max(c(nepoch,niter))))),' - Iteration ',formatC(ni,width=1+ceiling(log10(max(c(nepoch,niter))))),' - ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4),collapse=' - ')))
#         }else{
#           print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Epoch ',formatC(ne,width=1+ceiling(log10(max(c(nepoch,niter))))),' - Iteration ',formatC(ni,width=1+ceiling(log10(max(c(nepoch,niter)))))))
#         }
#       }
#     }
#   }
#
#   return(W_curr)
#
# }

# version where only a number of rows of the covariance matrix is estimated

# robinv <- function(A){
#   tryCatch(
#     expr={Ai <- inv(A)
#           return(Ai)},
#     warning=function(w){Ai <- solve.qr(qr(A))
#                         Ai[is.na(Ai)] <- 0
#                         return(Ai)},
#     error=function(e){print(e)}
#   )
# }

# robinv <- function(A){
#
#   Ai <- matrix(0,nrow=nrow(A),ncol=ncol(A))
#   diag(Ai) <- 1
#
#   tryCatch(
#     expr={Ai <- inv(A)},
#     warning=function(w){Ai <- solve.qr(qr(A),tol=1e-3)
#     Ai[is.na(Ai)] <- 0},
#     error=function(e){print(e)}
#   )
#
#   return(Ai)
# }

cpca_pimvr <- function(X_real,X_imag,W_init,nepoch,niter,beta,bs,nb,bsf,onr=1e2,verbose=F,prtint=1,W_ref=W_init,r2=F,z=F,wrkrs=1){

  # estimation of cPCA loading matrix based on the artice
  #   De Sa, Christopher, et al. "Accelerated Stochastic Power Iteration." arXiv preprint arXiv:1707.02670 (2017).
  #
  # inputs
  #   X_real  - (bigmemory) data matrix containing the real part of the data
  #   X_imag  - (bigmemory) data matrix containing the imaginary part of the data
  #   W_init  - intialization of the loading matrix, initialization serves also to devise the number of components to estimate
  #   nepoch  - number of epoch - in each epoch there is a costly pass over all the data
  #   niter   - number of iterations per epoch
  #   beta    - momentum parameter, either scalar or vector of length nepoch or nepoch*niter
  #   bs      - block size, size of the minibatches of data for estimation of the covariance matrix in each iteration
  #   nb      - number of blocks, number of the minibatches of data used for the estimation of the covariance matrix in each iteration
  #   bsf     - block size in full data pass, size of the batches of data for the estimation of the covariance matrix on the whole data at the start of each epoch
  #   onr     - number of rows of covariance matrix to estimate to reduce the memory usage by the algorithm
  #   verbose - logical, printing of progress
  #   prtint  - interval at which the information will be printed
  #   W_ref   - reference solution to the W estimation, only relevant when verbose=T, then shows modulus of correlation between W and W_ref
  #   r2      - logical, estimation of variance explained, printed if verbose=T
  #   z       - logical, should scores be returned
  #   wrkrs   - number of parallel workers
  #
  # outputs
  #   W  - resulting loading matrix
  #   R2 - explained variance per component, returned when r2=T


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

  if(r2){
    totvar_orig <- 0
    for(clmn in 1:ncol(X_real)){
      helper <- complex(real=X_real[,clmn],imaginary=X_imag[,clmn])
      totvar_orig <- totvar_orig + var(helper)
    }
  }

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

    plan(list(sequential,tweak(multicore,workers=wrkrs,gc=T,earlySignal=T)))

    hlpr_w <- listenv()
    rb <- c(seq(from=1,to=p,by=onr),p+1)
    for(ri in 1:(length(rb)-1)){
      hlpr_w[[ri]] %<-% {
        A_sub <- matrix(0,nrow=rb[ri+1]-rb[ri],ncol=p)
        hlpr_b <- listenv()
        for(bi in 1:nbf){
          hlpr_b[[bi]] %<-% {
            x1 <- matrix(
                    complex(real=     as.vector(X_real[bf_list[[bi]],rb[ri]:(rb[ri+1]-1)]),
                            imaginary=as.vector(X_imag[bf_list[[bi]],rb[ri]:(rb[ri+1]-1)])),
                          nrow=bsf,ncol=rb[ri+1]-rb[ri])
            x2 <- matrix(
                    complex(real=     as.vector(X_real[bf_list[[bi]],]),
                            imaginary=as.vector(X_imag[bf_list[[bi]],])),
                          nrow=bsf,ncol=p)
            Conj(t(x1)) %*% x2
          }
        }
        A_sub <- Reduce('+',as.list(hlpr_b))
        A_sub <- A_sub/nb
        A_sub %*% W_curr
      }
    }

    W_anch <- do.call(rbind,as.list(hlpr_w))

    tryCatch(
      expr={R_prev_inv <- pinv(R_prev)},
      warning=function(w){
                print(w)
                print(paste0('A problem occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                return(W_curr)},
      error=function(e){
              print(e)
              print(paste0('An error occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
              return(W_curr)})

    W_anch <- W_anch - beta[(ne-1)*niter+1]*t(t(R_prev_inv) %*% t(W_prev))
    R_anch <- qr.R(qr(rbind(W_anch,W_prev)))


    tryCatch(
      expr={R_anch_inv <- pinv(R_anch)},
      warning=function(w){
                print(w)
                print(paste0('A problem occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                return(W_curr)},
      error=function(e){
              print(e)
              print(paste0('An error occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
              return(W_curr)})

    W_anch <- W_anch %*% R_anch_inv

    hlpr_w <- listenv()
    rb <- c(seq(from=1,to=p,by=onr),p+1)
    for(ri in 1:(length(rb)-1)){
      hlpr_w[[ri]] %<-% {
        A_sub <- matrix(0,nrow=rb[ri+1]-rb[ri]+1,ncol=p)
        hlpr_b <- listenv()
        for(bi in 1:nbf){
          hlpr_b[[bi]] %<-% {
            x1 <- matrix(
                    complex(real=     as.vector(X_real[bf_list[[bi]],rb[ri]:(rb[ri+1]-1)]),
                           imaginary=as.vector(X_imag[bf_list[[bi]],rb[ri]:(rb[ri+1]-1)])),
                    nrow=bsf,ncol=rb[ri+1]-rb[ri])
            x2 <- matrix(
                    complex(real=     as.vector(X_real[bf_list[[bi]],]),
                            imaginary=as.vector(X_imag[bf_list[[bi]],])),
                    nrow=bsf,ncol=p)
            Conj(t(x1)) %*% x2
          }
        }
        A_sub <- Reduce('+',as.list(hlpr_b))
        A_sub <- A_sub/nb
        A_sub %*% W_anch
      }
    }
    V_anch <- do.call(rbind,as.list(hlpr_w))

    for(ni in 1:niter){

      b <- sample(x=1:n,size=bs*nb)
      b_list <- list()
      for(bi in 1:nb){
        b_list[[bi]] <- b[((bi-1)*bs+1):(bi*bs)]
      }

      W_prev_help <- W_curr

      alpha <- Conj(t(W_curr)) %*% W_anch

      hlpr_w <- listenv()
      rb <- c(seq(from=1,to=p,by=onr),p+1)
      for(ri in 1:(length(rb)-1)){
        hlpr_w[[ri]] %<-% {
          A_sub <- matrix(0,nrow=rb[ri+1]-rb[ri]+1,ncol=p)
          hlpr_b <- listenv()
          for(bi in 1:nb){
            hlpr_b[[bi]] %<-% {
              x1 <- matrix(
                      complex(real=     as.vector(X_real[b_list[[bi]],rb[ri]:(rb[ri+1]-1)]),
                              imaginary=as.vector(X_imag[b_list[[bi]],rb[ri]:(rb[ri+1]-1)])),
                      nrow=bs,ncol=rb[ri+1]-rb[ri])
              x2 <- matrix(
                      complex(real=     as.vector(X_real[b_list[[bi]],]),
                              imaginary=as.vector(X_imag[b_list[[bi]],])),
                      nrow=bs,ncol=p)
              Conj(t(x1)) %*% x2
            }
          }
          A_sub <- Reduce('+',as.list(hlpr_b))
          A_sub <- A_sub/nb
          A_sub %*% (W_curr - W_anch %*% alpha)
        }
      }
      W_curr <- do.call(rbind,as.list(hlpr_w))


      tryCatch(
        expr={R_prev_inv <- pinv(R_prev)},
        warning=function(w){
                  print(w)
                  print(paste0('A problem occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                  return(W_curr)},
        error=function(e){
                print(e)
                print(paste0('An error occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                return(W_curr)})

      W_curr <- W_curr + V_anch %*% alpha - beta[(ne-1)*niter+ni]*Conj(t(Conj(t(R_prev_inv)) %*% Conj(t(W_prev))))
      R_curr <- qr.R(qr(rbind(W_curr,W_prev_help)))

      tryCatch(
        expr={R_curr_inv <- pinv(R_curr)},
        warning=function(w){
                  print(w)
                  print(paste0('A problem occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                  return(W_curr)},
        error=function(e){
                print(e)
                print(paste0('An error occured, returning the current estimate at epoch ',ne,' and iteration ',ni,'.'))
                return(W_curr)})

      W_curr <- W_curr %*% R_curr_inv

      R_prev <- R_curr
      W_prev <- W_prev_help

      if(r2 | z){

        var_cpca <- rep(NA,times=ncol(W_curr))
        Z <- matrix(NA,nrow=nrow(X_real),ncol=ncol(W_curr))
        for(pc in 1:ncol(W_curr)){
          y <- listenv()
          for(no in 1:nrow(X_real)){y[[no]] %<-% {sum(complex(real=as.vector(X_real[no,]),imaginary=as.vector(X_imag[no,])) * W[,pc])}}
          Z[,pc] <- do.call(c,as.list(y))
        }
        explvari <- apply(Z,2,var)/totvar_orig
      }

      if(verbose & ((((ni - 1) %% prtint)) == 0)){
        if(!is.null(W_ref) | r2){cat(paste0('\n',paste0(rep('#',floor(getOption('width')/2)),collapse=' ')))}

        cat(paste0('\n',format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),' - Epoch ',formatC(ne,width=1+ceiling(log10(max(c(nepoch,niter))))),' - Iteration ',formatC(ni,width=1+ceiling(log10(max(c(nepoch,niter)))))))
        if(!is.null(W_ref)){
          mccs <- rep(NA,times=ncol(W_ref))
          for(u in 1:ncol(W_ref)){mccs[u] <- Mod(cor(x=W_ref[,u],y=W_curr[,u]))}
          cat(paste0('\n  CC: ',paste0('cPC',1:ncol(W_ref),': ',formatC(mccs,format='f',digits=4),collapse=' - ')))
        }
        if(r2){
          cat(paste0('\n  R2: ',paste0('cPC',1:length(explvari),': ',formatC(explvari,format='f',digits=4),collapse=' - ')))
        }
      }


    }
  }

  plan(sequential)

  if(r2 | z){
    if(z){
      return(list(W=W_curr,Z=Z,R2=explvari))
    }
    else{
      return(list(W=W_curr,R2=explvari))
    }
  }else{
    return(W_curr)
  }
}
