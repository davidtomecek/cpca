rm(list=ls(all.names=T))
graphics.off()

library(complexlm)
#
# get_explained_variance <- function(cpca,Data){
#
#   # computing individual explained variance from a group cPCA estimate
#   #
#   # inputs
#   # cpca ... estimated group cPCA
#   # Data ... list of individual data sets, the data sets has to be matched
#   #
#   # outputs
#   # explvari ... matrix of explained variance per component and data set
#
#   library(complexlm)
#
#   # extract individual loadings - individual waveforms corresponding to cPCA components
#   Loads_ind <- list()
#   for(j in 1:length(Data)){
#     m <- ncol(Data[[j]])
#     n <- sum(sapply(Data[1:(j-1)],ncol))
#     Loads_ind[[j]] <- cpca$rotation[(1+n):(m+n),]
#   }
#
#   # rescaling the indivdual loadings so that it's a unitary matrix
#   Loads_ind_resc <- list()
#   for(j in 1:length(Loads_ind)){
#     Loads_ind_resc[[j]] <- Loads_ind[[j]]
#     for(k in 1:ncol(Loads_ind_resc[[j]])){
#         Loads_ind_resc[[j]][,k] <- Loads_ind_resc[[j]][,k]/sqrt(sum(Mod(Loads_ind_resc[[j]][,k])^2))
#     }
#   }
#
#   # compute individual scores with rescaled loadings - individual maps of components on the brain surface
#   Scores_ind_resc <- list()
#   for(j in 1:length(Sets[[set]])){
#     k <- which(desc$dataset == Sets[[set]][j])
#     Scores_ind_resc[[j]] <- Data_comb[[k]] %*% Loads_ind_resc[[j]]
#   }
#
#   # explained variance from the group model on the level of each individual dataset
#   explvari <- matrix(NA,nrow=length(Data),ncol=ncol(Loads_ind_resc[[1]]))
#   for(j in 1:length(Data)){
#     var_data <- apply(Data[[j]],2,var)
#     var_scores <- apply(Scores_ind_resc[[j]],2,var)
#     explvari[j,] <- var_scores/sum(var_data)
#   }
#
#   return(explvari)
#
# }

# dec <- 'fsaverage5'
dec <- 'fsaverage6'

smooth <- 'sm5'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'

# # list files with results
# files_cpca <- list.files(path=path_work,pattern=paste0('workspace_',dec,'_cpca_pimvr_[0-9]{8}\\.RData'))
# or just load the specific results file
# load(paste0(path_work,'workspace_fsaverage5_cpca_pimvr_20251001.RData'))

# trying the intermediate results on the fsaverage6 data
# load(paste0(path_work,'workspace_fsaverage6_cpca_pimvr_partial.RData'))
load(paste0(path_work,'workspace_fsaverage6_sm5_cpca_pimvr_partial.RData'))

n_pc <- ncol(cPCA_list[[1]]$W)


ExplVari <- list()

for(set in 1:length(Sets)){

  explvari <- matrix(NA,nrow=length(Sets[[set]]),ncol=n_pc)
  id <- substr(Sets[[set]],5,10)
  group <- substr(Sets[[set]],5,5)
  visit <- substr(Sets[[set]],26,26)

#   # select the result with the greatest overall explained variance
#   R2 <- do.call(rbind,lapply(cPCA_list,function(x) x$R2))
#   R2_ord <- apply(R2,2,order,decreasing=T)
#
#   ind_best <- table(R2_ord[1,])
#   ind_best <- as.integer(names(ind_best)[which.max(ind_best)])
#
#   W_glob <- cPCA_list[[ind_best]]$W

  W_glob <- cPCA_list[[1]]$W

  # then go through the preprocessed files according to the files indicated in the Sets variable
  for(fold in 1:length(Sets[[set]])){

    load(paste0(path_work,Sets[[set]][fold],'/',dec,'_',smooth,'_group_cpca_preprocessed_data.RData'))

    data_preproc <- t(data_preproc[-coni,])

    scores_indi <- data_preproc %*% W_glob

    var_data <- apply(data_preproc,2,var)
    var_scores <- apply(scores_indi,2,var)
    explvari[fold,] <- var_scores/sum(var_data)

  }

  names(explvari) <- paste0('cPC',1:n_pc)
  explvari <- as.data.frame(explvari)
  explvari$id <- id
  explvari$group <- group
  explvari$visit <- visit

  ExplVari[[set]] <- explvari

}

# save(ExplVari,file=paste0(path_work,'workspace_fsaverage5_cpca_pimvr_20251001_explvari.RData'))
# save(ExplVari,file=paste0(path_work,'workspace_',dec,'_cpca_pimvr_partial_explvari.RData'))
save(ExplVari,file=paste0(path_work,'workspace_',dec,'_',smooth,'_cpca_pimvr_partial_explvari.RData'))

# ExplVari <- list()
# for(set in 1:length(Sets)){
#   ExplVari[[set]] <- get_explained_variance(cpca=CPCA[[set]],
#                                             Data=Data_comb[match(Sets[[set]],desc$dataset)])
# }
#
# save(list=c('desc','Sets','ExplVari'),file=paste0(path_work,dec,'_optnorm_low_explained_variance_individual_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
#
#
# # computing indiividual explained variance for the data normalized in the high frequency range
# load(paste0(path_work,tail(files_cpca_high,n=1)))
# load(paste0(path_work,tail(files_data_high,n=1)))
#
# ExplVari <- list()
# for(set in 1:length(Sets)){
#   ExplVari[[set]] <- get_explained_variance(cpca=CPCA[[set]],
#                                             Data=Data_comb[match(Sets[[set]],desc$dataset)])
# }
#
# save(list=c('desc','Sets','ExplVari'),file=paste0(path_work,dec,'_optnorm_high_explained_variance_individual_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
