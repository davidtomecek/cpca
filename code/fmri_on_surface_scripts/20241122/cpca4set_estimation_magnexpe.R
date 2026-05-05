rm(list=ls())
graphics.off()

library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)

library(freesurferformats)
library(gsignal)

library(rsvd)

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
path_save <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

dec <- 'fsaverage6'
path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')
path_labe <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/label/')

load(paste0(path_work,'workspace_preprocessed_data_',dec,'.RData'))


n_pc <- 5

Sets <- list(gcv1=desc$dataset[desc$group == 'c' & desc$visit == 1],
             gcv2=desc$dataset[desc$group == 'c' & desc$visit == 2],
             gpv1=desc$dataset[desc$group == 'p' & desc$visit == 1],
             gpv2=desc$dataset[desc$group == 'p' & desc$visit == 2],
             gc  =desc$dataset[desc$group == 'c'],
             gp  =desc$dataset[desc$group == 'p'],
             all =desc$dataset)

expe_par_coef <- c(2,5,10)
expe_par_n <- c(1,2,3)
expe_par_ntry <- 10

CPCA_magnexpe <- list()

for(set in 1){ # :length(Sets)){

  print(paste0('##### ',names(Sets)[set],' #####'))
  print(Sets[[set]])

  helper <- list()

  for(expe in 1:expe_par_ntry){

    # choosing random parameters of the experiment
    expe_n <- sample(x=expe_par_n,size=1)
    expe_ind <- sample(1:length(Sets[[set]]),size=expe_n)
    expe_coef <- sample(x=expe_par_coef,size=expe_n)

    # changing the magnitude of one or more datasets
    Data_expe <- Data_comb[which(desc$dataset %in% Sets[[set]])]
    for(dataset in 1:expe_n){Data_expe[[expe_ind[dataset]]] <- Data_expe[[expe_ind[dataset]]]*expe_coef[dataset]}
    data_expe <- do.call(cbind,Data_expe)

    print(paste0(Sys.time(),': Estimating cPCA - experiment ',expe))

    expe_cpca <- rpca(data_expe,k=n_pc,center=F,scale=F)

    helper[[expe]] <- list(expe_n=expe_n,
                           expe_ind=expe_ind,
                           expe_coef=expe_coef,
                           expe_cpca=expe_cpca)
  }

  CPCA_magnexpe[[set]] <- helper

  print(paste0(Sys.time(),': Finished set ',set,'.'))
}

print(paste0(Sys.time(),': Everything finished! Saving the results.'))
save('CPCA_magnexpe',file=paste0(path_work,'workspace_cpca4set_magnexpe_',dec,'_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
