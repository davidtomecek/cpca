# TODO

rm(list=ls())
graphics.off()

options(future.rng.onMisuse='ignore')

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
library(future)
library(listenv)

path_utils <- '/home/anyzjiri/NUDZ/fmri_on_surface/scripts/20250606/'
source(paste0(path_utils,'utils_pimvr.R'))

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

# how many times try the random initialization + cpca estimation
ntry <- 3

# parameters for the vrpim method
nepoch <- 4
niter <- 16
# beta <- 10^seq(from=-4,to=-10,length.out=niter)
beta <- 10^seq(from=-4,to=-10,length.out=nepoch)
bs <- 128 # block size
nb <- 32 # number of blocks
bsd <- 200 # block size per data set - for the anchor estimate
onr <- 250 # # size block of variables used in the computation of the covariance matrix in the block manner

Sets <- list()
Sets[[1]] <- fold_data[grepl(pattern='^ESO_C',x=fold_data) & grepl(pattern='1$',x=fold_data)]
# Sets[[1]] <- fold_data[grepl(pattern='^ESO_C',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]
# Sets[[1]] <- fold_data[grepl(pattern='^ESO_P',x=fold_data) & (grepl(pattern='1$',x=fold_data) | grepl(pattern='2$',x=fold_data))]
# Sets[[1]] <- fold_data[!(grepl(pattern='3392',x=fold_data) | grepl(pattern='3247',x=fold_data))]

Res <- list()

for(set in 1:length(Sets)){

  # (1) prepare the data

  Sizes <- list()
  Conis <- list()


  for(fold in 1:length(Sets[[set]])){

    if(!dir.exists(paste0(path_work,Sets[[set]][fold]))){dir.create(paste0(path_work,Sets[[set]][fold]),recursive=T)}

    if(file.exists(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData')) & !force_preproc){

      print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Using already preprocessed data for ',Sets[[set]][fold]))

      load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
      Sizes[[fold]] <- dim(data_preproc)
      Conis[[fold]] <- which(apply(data_preproc,1,function(x) sum((x - mean(x))^2)) == 0)
      rm('data_preproc')
      gc()

    }else{

      plan(multicore,workers=16)

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

      # checking for constant variables in each data set
      Conis[[fold]] <- which(apply(data_preproc,1,function(x) sum((x - mean(x))^2)) == 0)

      # saving the preprocessed data
      save(data_preproc,file=paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))

      # deleting the preprocessed data
      rm('data_preproc')
      gc()

      plan(sequential)

    }
  }


  # (2) construct the bigmemory matrices for real and imaginary parts

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Creating the bigmemory data matrices'))

  # size of the big.matrix
  data_size <- c(sum(sapply(Sizes,function(x) x[2])),Sizes[[1]][1])

  # indices of constant variables in every data set
  coni <- Reduce(intersect,Conis)
  data_size[2] <- data_size[2] - length(coni)

  if(length(list.files(path_work,pattern='^data_comb_')) > 0){unlink(paste0(path_work,list.files(path_work,pattern='^data_comb_')))}

  data_comb_real <- filebacked.big.matrix(nrow=data_size[1],ncol=data_size[2],backingpath=path_work,backingfile='data_comb_real',descriptorfile='data_comb_real.desc')
  data_comb_imag <- filebacked.big.matrix(nrow=data_size[1],ncol=data_size[2],backingpath=path_work,backingfile='data_comb_imag',descriptorfile='data_comb_imag.desc')

  curr_row <- 1
  for(fold in 1:length(Sets[[set]])){

    load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))

    dimnames(data_preproc) <- NULL

    data_comb_real[curr_row:(curr_row+Sizes[[set]][2]-1),] <- Re(t(data_preproc)[,-coni])
    data_comb_imag[curr_row:(curr_row+Sizes[[set]][2]-1),] <- Im(t(data_preproc)[,-coni])

    curr_row <- curr_row + Sizes[[set]][2]
  }

  # (3) centering and scaling the data

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Centering and scaling the data'))

  # centering and scaling the data
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
#   for(clmn in 1:ncol(data_comb_real)){
#     helper <- complex(real=data_comb_real[,clmn],imaginary=data_comb_imag[,clmn])
#
#     mu <- mean(helper)
#     sigma2 <- var(helper)
#     helper <- (helper - mu)/sqrt(sigma2)
#
#     data_comb_real[,clmn] <- Re(helper)
#     data_comb_imag[,clmn] <- Im(helper)
#   }

  # implicit parallelism with BLAS
  blas_set_num_threads(16)
  omp_set_num_threads(16)

 # (X) estimate the reference solution using rpca from the rsvd package

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating reference group cPCA using rpca from package rsvd'))
  data_comb <- matrix(complex(real=as.matrix(data_comb_real),imaginary=as.matrix(data_comb_imag)),nrow=data_size[1],ncol=data_size[2])
  cpca_opt <- rpca(A=data_comb,k=n_pc,center=F,scale=F)
  rm('data_comb')
  gc()

  # (4) VRPIM method in group setting

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating group cPCA using the PIMVR method'))

  cPCA_list <- list()

  for(nt in 1:ntry){

    # (4.1) initialization from a small sample of the data

    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating initialization on a subset of data'))

    ind_init <- sample(1:nrow(data_comb_real),size=size_init)
    data_init <- matrix(complex(real=as.matrix(data_comb_real[ind_init,]),imaginary=as.matrix(data_comb_imag[ind_init,])),nrow=size_init,ncol=ncol(data_comb_real))

    cpca_init <- rpca(A=data_init,k=n_pc,center=F,scale=F)
    W <- cpca_init$rotation

    rm('ind_init','data_init')
    gc()

    # (4.2) estimating group cpca by the VRPIM method

    blas_set_num_threads(16)
    omp_set_num_threads(16)

    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating of the cPCA using the PIMVR method'))


#     W <- cpca_pim(X_real=data_comb_real,
#                   X_imag=data_comb_imag,
#                   W_init=W,
#                   niter=niter,
#                   beta=beta,
#                   bs=bs,
#                   nb=nb,
#                   onr=onr,
#                   verbose=T,
#                   prtint=1,
#                   W_ref=cpca_opt$rotation,
#                   wrkrs=4)

    bsf <- length(Sets[[set]])*bsd # block size for the anchor estimate
    if(bsf > data_size[1]){stop('The block size for the anchor covariance matrix estimate is larger than the number of observations.')}

    W <- cpca_pimvr(X_real=data_comb_real,
                    X_imag=data_comb_imag,
                    W_init=W,
                    nepoch=nepoch,
                    niter=niter,
                    beta=beta,
                    bs=bs,
                    nb=nb,
                    bsf=bsf,
                    onr=onr,
                    verbose=T,
                    prtint=1,
                    W_ref=cpca_opt$rotation,
                    wrkrs=16)

    cPCA_list[[nt]] <- list(W=W)

    # (5) estimating the explained variance

    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Estimating the variance explained'))

    var_orig <- rep(NA,times=ncol(data_comb_real))
    for(clmn in 1:ncol(data_comb_real)){
      helper <- complex(real=data_comb_real[,clmn],imaginary=data_comb_imag[,clmn])
      var_orig[clmn] <- var(helper)
    }

    var_cpca <- rep(NA,times=n_pc)
    for(pc in 1:ncol(W)){

      plan(multicore,workers=16)
      z <- foreach(no=1:nrow(data_comb_real),.combine='c') %dofuture% {
             sum(complex(real=as.vector(data_comb_real[no,]),imaginary=as.vector(data_comb_imag[no,])) * W[,pc])
           }
      plan(sequential)
      var_cpca[pc] <- var(z)
    }

    R2 <- var_cpca/sum(var_orig)

    print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Variance explained in the ',nt,' run'))
    print(formatC(R2,format='f',digits=3,width=3))

    cPCA_list[[nt]]$R2 <- R2

  }

  print(paste0(format(Sys.time(),format='%Y-%m-%d %H:%M:%S'),': Summary of the variance explained in the ',ntry,' runs'))

  R2 <- do.call(rbind,lapply(cPCA_list,function(x) x$R2))
  print(R2)
  print(apply(R2,2,order))

  Res[[set]] <- list(coni=coni,
                     n_pc=n_pc,
                     size_init=size_init,
                     ntry=ntry,
#                      nepoch=nepoch,
                     niter=niter,
                     beta=beta,
                     bs=bs,
                     nb=nb,
#                      bsd=bsd,
                     onr=onr,
                     Sizes=Sizes,
                     cPCA_list=cPCA_list)

}

save(list=c('Res','Sets','dec','coni','n_pc','size_init','ntry','niter','beta','bs','nb','onr'), # 'nepoch','bsd',
     file=paste0(path_save,'workspace_',dec,'_cpca_pimvr_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))

blas_set_num_threads(1)
omp_set_num_threads(1)
