rm(list=ls())
graphics.off()

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)

library(rsvd)

# library(doParallel)
# compclust <- makeForkCluster(20)
# registerDoParallel(compclust)

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
    m <- lm(data_left[row,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    x <- resid(m)
    return(x)
  }

  data_right <- foreach(row=1:nrow(data_right),.combine='rbind') %dofuture% {
    m <- lm(data_right[row,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
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

path_data <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/rest/'
fold_data <- list.dirs(path_data,full.names=F,recursive=F)

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

n_pc <- 5

visit <- as.integer(substr(fold_data,nchar(fold_data[1]),nchar(fold_data[1])))
group <- tolower(substr(fold_data,5,5))

Sets <- list(gcv1=fold_data[group == 'c' & visit == 1],
             gcv2=fold_data[group == 'c' & visit == 2],
             gpv1=fold_data[group == 'p' & visit == 1],
             gpv2=fold_data[group == 'p' & visit == 2],
             gc  =fold_data[group == 'c'],
             gp  =fold_data[group == 'p'],
             all =fold_data)

# plan(sequential)
plan(multisession,workers=16)

Data_comb <- foreach(fold=1:length(fold_data)) %do% {

  data <- preprocess_dataset(path_data_left=paste0(path_data,fold_data[fold],'/rest/001/fmcpr.sm0.',dec,'.lh.nii.gz'),
                             path_data_right=paste0(path_data,fold_data[fold],'/rest/001/fmcpr.sm0.',dec,'.rh.nii.gz'),
                             path_mcp=paste0(path_data,fold_data[fold],'/rest/001/fmcpr.mcdat'),
                             subset_left=veri_left,
                             subset_right=veri_right,
                             fs_orig=1,
                             fs_up=2)
  return(data)
}

plan(sequential)

CPCA <- list()

for(set in 1:length(Sets)){

  print(paste0('##### ',names(Sets)[set],' #####'))
  print(Sets[[set]])
  print(paste0(Sys.time(),': Preprocessing data'))

  data_comb <- do.call(cbind,Data_comb[which(fold_data %in% Sets[[set]])])

  print(paste0(Sys.time(),': Estimating cPCA'))

  CPCA[[set]] <- rpca(data_comb,k=n_pc,center=F,scale=F)

}

desc <- data.frame(dataset=fold_data,
                   group=group,
                   visit=visit)

print(paste0(Sys.time(),': Finished! Saving the results.'))
save(list=c('hemi_left','hemi_right','veri_left','veri_right','Sets','desc','Data_comb'),
     file=paste0(path_work,'workspace_cpca4set_optnorm_low_data_',dec,'_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
save(list=c('hemi_left','hemi_right','veri_left','veri_right','Sets','desc','CPCA'),
     file=paste0(path_work,'workspace_cpca4set_optnorm_low_cpca_',dec,'_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
