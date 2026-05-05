rm(list=ls(all.names=T))
graphics.off()

library(complexlm)

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
files <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9]{8}\\.RData$'))

load(paste0(path_work,tail(files,n=1)))


ExplVari <- list()

for(set in 1:length(Sets)){
# for(set in 1){

  # extract individual loadings - individual waveforms corresponding to cPCA components
  Loads_ind <- list()
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    if(j == 1){
      m <- ncol(Data_comb[[k]])
      Loads_ind[[j]] <- CPCA[[set]]$rotation[1:m,]
    }else{
      m <- ncol(Data_comb[[k]])
      n <- sum(sapply(Data_comb[which(desc$dataset[1:(j-1)] %in% desc$dataset)],ncol))
      Loads_ind[[j]] <- CPCA[[set]]$rotation[(1+n):(m+n),]
    }
  }

  # rescaling the indivdual loadings so that it's a unitary matrix
  Loads_ind_resc <- list()
  for(j in 1:length(Loads_ind)){
    Loads_ind_resc[[j]] <- Loads_ind[[j]]
    for(k in 1:ncol(Loads_ind_resc[[j]])){
        Loads_ind_resc[[j]][,k] <- Loads_ind_resc[[j]][,k]/sqrt(sum(Mod(Loads_ind_resc[[j]][,k])^2))
    }
  }

  # compute individual scores - individual maps of components on the brain surface
  Scores_ind <- list()
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    Scores_ind[[j]] <- Data_comb[[k]] %*% Loads_ind[[j]]
  }

  # and also the individual scores with rescaled loadings
  Scores_ind_resc <- list()
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    Scores_ind_resc[[j]] <- Data_comb[[k]] %*% Loads_ind_resc[[j]]
  }

  # rearrange individual scores by cPC to get the scores in the group form
  Scores_group <- list()
  for(j in 1:ncol(Scores_ind[[1]])){
    Scores_group[[j]] <- list(pc=j,
                              scores=do.call(rbind,lapply(Scores_ind,function(x) t(x[,j]))))
  }


  # explained variance from the group model on the level of each individual dataset
  explvari <- matrix(NA,nrow=length(Sets[[set]]),ncol=ncol(Loads_ind[[1]]))
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    var_data <- apply(Data_comb[[k]],2,var)
    var_scores <- apply(Scores_ind_resc[[j]],2,var)
    explvari[j,] <- var_scores/sum(var_data)
  }

  ExplVari[[set]] <- explvari

}

save(list=c('desc','Sets','ExplVari'),file=paste0(path_work,dec,'_explained_variance_individual_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
