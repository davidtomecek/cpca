rm(list=ls(all.names=T))
graphics.off()

get_explained_variance <- function(cpca,Data){

  # computing individual explained variance from a group cPCA estimate
  #
  # inputs
  # cpca ... estimated group cPCA
  # Data ... list of individual data sets, the data sets has to be matched
  #
  # outputs
  # explvari ... matrix of explained variance per component and data set

  library(complexlm)

  # extract individual loadings - individual waveforms corresponding to cPCA components
  Loads_ind <- list()
  for(j in 1:length(Data)){
    m <- ncol(Data[[j]])
    n <- sum(sapply(Data[1:(j-1)],ncol))
    Loads_ind[[j]] <- cpca$rotation[(1+n):(m+n),]
  }

  # rescaling the indivdual loadings so that it's a unitary matrix
  Loads_ind_resc <- list()
  for(j in 1:length(Loads_ind)){
    Loads_ind_resc[[j]] <- Loads_ind[[j]]
    for(k in 1:ncol(Loads_ind_resc[[j]])){
        Loads_ind_resc[[j]][,k] <- Loads_ind_resc[[j]][,k]/sqrt(sum(Mod(Loads_ind_resc[[j]][,k])^2))
    }
  }

  # compute individual scores with rescaled loadings - individual maps of components on the brain surface
  Scores_ind_resc <- list()
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    Scores_ind_resc[[j]] <- Data_comb[[k]] %*% Loads_ind_resc[[j]]
  }

  # explained variance from the group model on the level of each individual dataset
  explvari <- matrix(NA,nrow=length(Data),ncol=ncol(Loads_ind_resc[[1]]))
  for(j in 1:length(Data)){
    var_data <- apply(Data[[j]],2,var)
    var_scores <- apply(Scores_ind_resc[[j]],2,var)
    explvari[j,] <- var_scores/sum(var_data)
  }

  return(explvari)

}

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

# files for data normalized in f <= 0.15 Hz part of the spectrum
files_cpca_low <- list.files(path_work,pattern=paste0('^workspace_cpca4set_optnorm_low_cpca_',dec,'_[0-9]{8}\\.RData$'))
files_data_low <- list.files(path_work,pattern=paste0('^workspace_cpca4set_optnorm_low_data_',dec,'_[0-9]{8}\\.RData$'))

# files for data normalized in f > 0.15 Hz part of the spectrum
files_cpca_high <- list.files(path_work,pattern=paste0('^workspace_cpca4set_optnorm_high_cpca_',dec,'_[0-9]{8}\\.RData$'))
files_data_high <- list.files(path_work,pattern=paste0('^workspace_cpca4set_optnorm_high_data_',dec,'_[0-9]{8}\\.RData$'))


# computing indiividual explained variance for the data normalized in the low frequency range
load(paste0(path_work,tail(files_cpca_low,n=1)))
load(paste0(path_work,tail(files_data_low,n=1)))

ExplVari <- list()
for(set in 1:length(Sets)){
  ExplVari[[set]] <- get_explained_variance(cpca=CPCA[[set]],
                                            Data=Data_comb[match(Sets[[set]],desc$dataset)])
}

save(list=c('desc','Sets','ExplVari'),file=paste0(path_work,dec,'_optnorm_low_explained_variance_individual_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))


# computing indiividual explained variance for the data normalized in the high frequency range
load(paste0(path_work,tail(files_cpca_high,n=1)))
load(paste0(path_work,tail(files_data_high,n=1)))

ExplVari <- list()
for(set in 1:length(Sets)){
  ExplVari[[set]] <- get_explained_variance(cpca=CPCA[[set]],
                                            Data=Data_comb[match(Sets[[set]],desc$dataset)])
}

save(list=c('desc','Sets','ExplVari'),file=paste0(path_work,dec,'_optnorm_high_explained_variance_individual_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
