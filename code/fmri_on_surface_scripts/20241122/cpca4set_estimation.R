rm(list=ls())
graphics.off()

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)

library(doMC)
registerDoMC(10)

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

  data_left <- foreach(j=1:nrow(data_left),.combine='rbind') %dopar% {
    m <- lm(data_left[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
    x <- resid(m)
    x <- approx(t_orig,x,t_up)$y
    x <- get_raut(x=x,fs=fs_up,fsdec=fs_up)
    return(x)
  }

  data_right <- foreach(j=1:nrow(data_right),.combine='rbind') %dopar% {
    m <- lm(data_right[j,] ~ crpar[,1] + crpar[,2] + crpar[,3] + crpar[,4] + crpar[,5] + crpar[,6])
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

# number of principal components to retain
n_pc <- 5

visit <- as.integer(substr(fold_data,nchar(fold_data[1]),nchar(fold_data[1])))
group <- tolower(substr(fold_data,5,5))

Sets <- list(gcv1=fold_data[group == 'c' & visit == 1],
             gcv2=fold_data[group == 'c' & visit == 2],
             gpv1=fold_data[group == 'p' & visit == 1],
             gpv2=fold_data[group == 'p' & visit == 2])

CPCA <- list()

for(set in 1:length(Sets)){

  print(paste0('##### ',names(Sets)[set],' #####'))
  print(Sets[[set]])
  print(paste0(Sys.time(),': Preprocessing data'))

  data_comb <- foreach(fold=1:length(Sets[[set]]),.combine='cbind') %dopar% {
  
    data <- preprocess_dataset(path_data_left=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.sm0.',dec,'.lh.nii.gz'),
                               path_data_right=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.sm0.',dec,'.rh.nii.gz'),
                               path_mcp=paste0(path_data,Sets[[set]][fold],'/rest/001/fmcpr.mcdat'),
                               subset_left=veri_left,
                               subset_right=veri_right,
                               fs_orig=1,
                               fs_up=3)

    return(data)

  }

  print(paste0(Sys.time(),': Estimating cPCA'))

  CPCA[[set]] <- prcomp(data_comb,center=F,scale=F)
  CPCA[[set]]$rotation <- CPCA[[set]]$rotation[,1:n_pc]
  CPCA[[set]]$center <- NULL
  CPCA[[set]]$scale <- NULL
  CPCA[[set]]$x <- CPCA[[set]]$x[,1:n_pc]

}

print(paste0(Sys.time(),': Finished! Saving the results.'))
save(list=c('hemi_left','hemi_right','veri_left','veri_right','Sets','CPCA'),file=paste0(path_work,'workspace_cpca4set_',dec,'_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
