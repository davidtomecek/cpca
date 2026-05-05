rm(list=ls(all.names=T))
graphics.off()

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
files <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9]{8}\\.RData$'))

load(paste0(path_work,tail(files,n=1)))

# setting up what to test in the different group analyses

# Sets <- list(gcv1=fold_data[group == 'c' & visit == 1],
#              gcv2=fold_data[group == 'c' & visit == 2],
#              gpv1=fold_data[group == 'p' & visit == 1],
#              gpv2=fold_data[group == 'p' & visit == 2],
#              gc  =fold_data[group == 'c'],
#              gp  =fold_data[group == 'p'],
#              all =fold_data)

Tests <- list()
Tests[[1]] <- NULL # healthy controls at V1 -> nothing to test
Tests[[2]] <- NULL # healthy controls at V2 -> nothing to test
Tests[[3]] <- NULL # patients at V1 -> nothing to test
Tests[[4]] <- NULL # patients at V2 -> nothing to test
# healthy controls at V1 and V2
setind <- which(Sets[[5]] %in% desc$dataset)
Tests[[5]] <- list(gcv2_gcv1=cbind(desc$visit[setind] == 2,desc$visit[setind] == 1))
# patients at V1 and V2
setind <- which(Sets[[6]] %in% desc$dataset)
Tests[[6]] <- list(gpv2_gpv1=cbind(desc$visit[setind] == 2,desc$visit[setind] == 1))
# all subjects at V1 and V2
setind <- which(Sets[[7]] %in% desc$dataset)
Tests[[7]] <- list(gcv2_gcv1=cbind(desc$visit[setind] == 2 & desc$group[setind]  == 'c',desc$visit[setind] == 1 & desc$group[setind]  == 'c'),
                   gpv2_gpv1=cbind(desc$visit[setind] == 2 & desc$group[setind]  == 'p',desc$visit[setind] == 1 & desc$group[setind]  == 'p'),
                   gpv1_gcv1=cbind(desc$visit[setind] == 1 & desc$group[setind]  == 'p',desc$visit[setind] == 1 & desc$group[setind]  == 'c'),
                   gpv2_gcv2=cbind(desc$visit[setind] == 2 & desc$group[setind]  == 'p',desc$visit[setind] == 2 & desc$group[setind]  == 'c'),
                   v2_v1    =cbind(desc$visit[setind] == 2,                             desc$visit[setind] == 1),
                   gp_gc    =cbind(                          desc$group[setind]  == 'p',                          desc$group[setind]  == 'p'))

for(set in 1:length(Sets)){
# for(set in 5){

  if(is.null(Tests[[set]])){next}

  # extract individual loadings
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

  # compute individual scores
  Scores_ind <- list()
  for(j in 1:length(Sets[[set]])){
    k <- which(desc$dataset == Sets[[set]][j])
    Scores_ind[[j]] <- Data_comb[[k]] %*% Loads_ind[[j]]
  }

  # rearrange individual scores by cPC
  Scores_group <- list()
  for(j in 1:ncol(Scores_ind[[1]])){
    Scores_group[[j]] <- list(pc=j,
                              scores=do.call(rbind,lapply(Scores_ind,function(x) t(x[,j]))))
  }

  # test the scores in various settings
  # starting with wilcoxon test
  P_values <- list()
  for(test in 1:length(Tests[[set]])){
    P_values[[test]] <- list()
    P_values[[test]]$test <- names(Tests[[set]])[test]
    for(pc in 1:length(Scores_group)){
      P_values[[test]][[paste0('p_value_pc',formatC(pc,width=3,flag=0),'_mod')]] <- apply(Scores_group[[pc]]$scores,2,function(x,ind){wilcox.test(Mod(x[ind[,1]]),Mod(x[ind[,2]]),exact=F)$p.value},ind=Tests[[set]][[test]])
      P_values[[test]][[paste0('p_value_pc',formatC(pc,width=3,flag=0),'_arg')]] <- apply(Scores_group[[pc]]$scores,2,function(x,ind){wilcox.test(Arg(x[ind[,1]]),Arg(x[ind[,2]]),exact=F)$p.value},ind=Tests[[set]][[test]])

      P_values[[test]][[paste0('mean_diff_pc',formatC(pc,width=3,flag=0),'_mod')]] <- apply(Scores_group[[pc]]$scores,2,function(x,ind){mean(Mod(x[ind[,1]])) - mean(Mod(x[ind[,2]]))},ind=Tests[[set]][[test]])
      P_values[[test]][[paste0('mean_diff_pc',formatC(pc,width=3,flag=0),'_arg')]] <- apply(Scores_group[[pc]]$scores,2,function(x,ind){mean(Arg(x[ind[,1]])) - mean(Arg(x[ind[,2]]))},ind=Tests[[set]][[test]])
    }
  }

  Test <- Tests[[set]]
  set <- names(Sets)[set]

  # saving the results
  save(list=c('hemi_left','hemi_right','veri_left','veri_right','set','Test','P_values'),
       file=paste0(path_work,'workspace_cpca4set_',dec,'_test_',set,'.RData'))

}
