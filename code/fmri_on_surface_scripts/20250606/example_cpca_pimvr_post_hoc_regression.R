rm(list=ls(all.names=T))
graphics.off()


visualize_on_surface <- function(x,surf,file,cmap='Zissou 1',ncols=256,range=NULL){

  cmap <- hcl.colors(n=ncols,palette=cmap)

  if(is.null(range)){
    cols <- cmap[round((ncols-1)*(x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))+1)]
    cols[is.na(x)] <- 'gray66'
  }else{
    cols <- cmap[round((ncols-1)*(x - range[1])/(range[2] - range[1])+1)]
    cols[is.na(x)] <- 'gray66'
  }

  suppressMessages({

    library(rgl)
    dev <- open3d()
    clear3d(type='lights')
    light3d(ambient='gray33',specular='gray33')
    shade3d(mesh3d(vertices=rbind(t(surf$vertices),rep(1,times=nrow(surf$vertices))),
                  triangles=t(surf$faces),
                  material=list(color=cols)),
              meshColor='vertices')
    widget <- rglwidget(width=1e3,height=1e3)
    htmlwidgets::saveWidget(widget=widget,file=file,selfcontained=T)
    close3d(dev=dev)
    rgl.quit()

  })

}

visualize_on_surface2 <- function(x,surf,cols,file){

  library(rgl)
  dev <- open3d()
  clear3d(type='lights')
  light3d(ambient='gray33',specular='gray33')
  shade3d(mesh3d(vertices=rbind(t(surf$vertices),rep(1,times=nrow(surf$vertices))),
                 triangles=t(surf$faces),
                 material=list(color=cols)),
            meshColor='vertices')
  widget <- rglwidget(width=1e3,height=1e3)
  htmlwidgets::saveWidget(widget=widget,file=file,selfcontained=T)
  close3d(dev=dev)
  rgl.quit()

}

pval_sig_cols <- function(x){
  cols <- rep('gray66',times=length(x))
  cols[x < 0.100] <- 'firebrick'
  cols[x < 0.050] <- 'orange'
  cols[x < 0.010] <- 'yellow'
  cols[x < 0.001] <- 'green'
  return(cols)
}

compute_cluster_size <- function(sigind,neigh,minsize=0,verbose=T){

  # function to compute cluster sizes
  #
  # input
  #   sigind  - indices of significant vertices
  #   neigh   - vertex neighborhood
  #   minsize - minimal size for a cluster
  #
  # output
  #   Clusts - list of clusters containing

  clust <- 1
  Clusts <- list()

  sigind_work <- sigind

  while(length(sigind_work) > 0){

    clustinds <- sigind_work[1]
    clustinds_open <- sigind_work[1]
    if(length(sigind_work) > 1){sigind_work <- sigind_work[2:length(sigind_work)]}else{break}

    while(length(clustinds_open) > 0){

      for(vert in 1:length(neigh[[clustinds_open[1]]])){

        matchind <- which(sigind_work %in% neigh[[clustinds_open[1]]])

        if(length(matchind) > 0){

          clustinds <- c(clustinds,sigind_work[matchind])
          clustinds_open <- c(clustinds_open,sigind_work[matchind])
          clustinds_open <- clustinds_open[2:length(clustinds_open)]

          sigind_work <- sigind_work[-matchind]

        }
      }

      if(length(clustinds_open) > 1){clustinds_open <- clustinds_open[2:length(clustinds_open)]}else{break}

    }

    Clusts[[clust]] <- list(indi=clustinds,
                            size=length(clustinds))
    clust <- clust + 1

  }

  for(clust in rev(which(sapply(Clusts,function(x) x$size) < minsize))){Clusts[[clust]] <- NULL}

  if(verbose){
    cat('Done, found ',length(Clusts),'clusters.\n')
  }

  return(Clusts)

}

wrap_t.test <- function(x,v) {
  p <- NA
  try({p <- t.test(x ~ v)$p.value},
      silent=T)
  return(p)
}

wrap_paired_t.test <- function(x,v) {
  p <- NA
#   try({p <- t.test(x ~ v,paired=T)$p.value}) # ,silent=T
  try({p <- t.test(x[v == levels(v)[1]],x[v == levels(v)[2]],paired=T)$p.value},
      silent=T)
  return(p)
}

wrap_wilcox.test <- function(x,v) {
  p <- NA
  try({p <- wilcox.test(x ~ v)$p.value},
      silent=T)
  return(p)
}

wrap_paired_wilcox.test <- function(x,v) {
  p <- NA
  try({p <- wilcox.test(x[v == levels(v)[1]],x[v == levels(v)[2]],paired=T)$p.value},
      silent=T)
  return(p)
}

permute_rep_meas <- function(lab,id){
  labspl <- split(lab,id)
  labspluni <- sapply(labspl,unique)
  labspluniperm <- sample(labspluni)
  labsplperm <- unsplit(mapply(function(x,y) rep(y,length(x)),labspl,labspluniperm,SIMPLIFY=F),id)
  return(labsplperm)
}

library(complexlm)
library(future)
library(listenv)

dec <- 'fsaverage6'

smooth <- 'sm5'

path_data <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'
# folds <- list.dirs(path_data,full.names=F,recursive=F)

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'

load(paste0(path_work,'workspace_',dec,'_',smooth,'_cpca_pimvr_partial.RData'))
# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_20251114.RData'))


setind <- 1
cpcaind <- 1

# plan(multicore,workers=16)
#
# print('Post-hoc regression')
#
# Betas <- listenv()
# for(fold in 1:length(Sets[[setind]])){
#
#   Betas[[fold]] %<-%{
#     load(paste0(path_data,Sets[[setind]][fold],'/',dec,'_',smooth,'_group_cpca_preprocessed_data.RData'))
#
#     data_preproc <- data_preproc[-coni,]
#
# #     Z <- Conj(t(data_preproc)) %*% Res[[setind]]$cPCA_list[[cpcaind]]$W
#     Z <- Conj(t(data_preproc)) %*% cPCA_list[[cpcaind]]$W
#
# #     B <- matrix(NA,nrow=nrow(Res[[setind]]$cPCA_list[[cpcaind]]$W),ncol=ncol(Res[[setind]]$cPCA_list[[cpcaind]]$W))
#     B <- matrix(NA,nrow=nrow(cPCA_list[[cpcaind]]$W),ncol=ncol(cPCA_list[[cpcaind]]$W))
#     for(pc in 1:ncol(Z)){
#       for(vert in 1:nrow(data_preproc)){
#         B[vert,pc] <- coef(lm(y ~ x,data=data.frame(y=data_preproc[vert,],x=Z[,pc])))[2]
#   #       B[vert,pc] <- coef(lm(y ~ x,data=as.data.frame(lapply(data.frame(y=data_preproc[vert,],x=Z[,pc]),scale,center=T,scale=T))))[2]
#
#       }
#     }
#
#     return(B)
#   }
# }
#
# Betas <- as.list(Betas)
#
# plan(sequential)
#
# save('Betas',file=paste0(path_work,'workspace_',dec,'_',smooth,'_cpca_pimvr_post_hoc_betas_20251208.RData'))

# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_post_hoc_betas_20251029.RData'))
# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_post_hoc_betas_20251030.RData'))

# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_post_hoc_betas_20251111.RData')) # regression coefficient
# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_post_hoc_betas_20251111.RData')) # correlation coefficient

# load(paste0(path_work,'workspace_',dec,'_cpca_pimvr_post_hoc_betas_20251118.RData'))

load(paste0(path_work,'workspace_fsaverage6_sm5_cpca_pimvr_post_hoc_betas_20251208.RData'))

# example - difference HC V1 + V2 x P V1 in cPC2
path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')
path_labe <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/label/')

path_rend <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/'
foldname_rend_ind <- 'cpa_global_pimvr_smooth_post_hoc_reg/individual'
# if(!dir.exists(paste0(path_rend,'/',foldname_rend_ind,'/'))){dir.create(paste0(path_rend,'/',foldname_rend_ind,'/'),recursive=T)}
unlink(paste0(path_rend,'/',foldname_rend_ind,'/'),recursive=T)
dir.create(paste0(path_rend,'/',foldname_rend_ind,'/'),recursive=T)

foldname_rend_comp <- 'cpa_global_pimvr_smooth_post_hoc_reg/comparison'
# if(!dir.exists(paste0(path_rend,'/',foldname_rend_comp,'/'))){dir.create(paste0(path_rend,'/',foldname_rend_comp,'/'),recursive=T)}
unlink(paste0(path_rend,'/',foldname_rend_comp,'/'),recursive=T)
dir.create(paste0(path_rend,'/',foldname_rend_comp,'/'),recursive=T)

library(freesurferformats)
inflated_left  <- read.fs.surface(paste0(path_surf,'lh.inflated'))
inflated_right <- read.fs.surface(paste0(path_surf,'rh.inflated'))

atlas_left  <- read.fs.annot(paste0(path_labe,'lh.aparc.annot'))
atlas_right <- read.fs.annot(paste0(path_labe,'rh.aparc.annot'))

hemi_left  <- 1:length(atlas_left$vertices)
hemi_right <- (1:length(atlas_left$vertices)) + length(atlas_left$vertices)

veri_left  <-  which(!(atlas_left$label_names %in% c('unknown','corpuscallosum')))
veri_right <- which(!(atlas_right$label_names %in% c('unknown','corpuscallosum')))

ind_left  <- 1:(length(veri_left))
ind_right <- (1+length(veri_left)):(length(veri_left) + length(veri_right))

rm('atlas_left','atlas_right')
gc()

# adjusting the sizes ...
n_pc <- ncol(Betas[[1]])
Betas <- lapply(Betas,
                function(x) {
                  y <- matrix(NA,nrow=Sizes[[1]][1],ncol=n_pc)
                  y[-coni,] <- x
                  return(y)
                })

# neighborhood lists
path_neigh <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/',dec,'/')

load(paste0(path_neigh,'neigh_left.RData'))
load(paste0(path_neigh,'neigh_right.RData'))


# visualizing individual maps based on post-hoc regression

n_pc <- 3

# for(set in seq(from=1,to=length(Sets[[1]]),by=20)){
for(set in sample(x=1:length(Sets[[1]]),size=3)){

  for(pc in 1:n_pc){

#     library(rgl)

    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- sqrt(Mod(Betas[[set]][ind_left,pc]))
    visualize_on_surface(x=helper,
                        surf=inflated_left,
                        file=paste0(path_rend,'/',foldname_rend_ind,'/',Sets[[1]][set],'_',dec,'_',smooth,'_cpca_global_pimvr_left_mod_cpc',pc,'.html'))

    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- sqrt(Mod(Betas[[set]][ind_right,pc]))
    visualize_on_surface(x=helper,
                        surf=inflated_right,
                        file=paste0(path_rend,'/',foldname_rend_ind,'/',Sets[[1]][set],'_',dec,'_',smooth,'_cpca_global_pimvr_right_mod_cpc',pc,'.html'))

    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- Arg(Betas[[set]][ind_left,pc])
    visualize_on_surface(x=helper,
                        surf=inflated_left,
                        file=paste0(path_rend,'/',foldname_rend_ind,'/',Sets[[1]][set],'_',dec,'_',smooth,'_cpca_global_pimvr_left_arg_cpc',pc,'.html'))

    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- Arg(Betas[[set]][ind_right,pc])
    visualize_on_surface(x=helper,
                        surf=inflated_right,
                        file=paste0(path_rend,'/',foldname_rend_ind,'/',Sets[[1]][set],'_',dec,'_',smooth,'_cpca_global_pimvr_right_arg_cpc',pc,'.html'))


#     rgl.quit()
  }



}

visit <- substr(Sets[[1]],nchar(Sets[[1]][1]),nchar(Sets[[1]][1]))
group <- substr(Sets[[1]],5,5)
ident <- tolower(substr(Sets[[1]],5,10))

Res <- list()

for(pc in 1:n_pc){

  betas_cv1v2_pv1 <- sqrt(
                        Mod(
                          rbind(
                            do.call(rbind,lapply(Betas[group == 'C'               ],function(x) x[,pc])),
                            do.call(rbind,lapply(Betas[group == 'P' & visit == '1'],function(x) x[,pc]))
                          )
                        )
                      )
  labs_cv1v2_pv1 <- factor(c(rep('c_v1v2',sum(group=='C')),rep('p_v1',sum(group=='P' & visit=='1'))),levels=c('c_v1v2','p_v1'))
  ids_cv1v2_pv1 <- c(ident[group == 'C'],ident[group == 'P' & visit == '1'])

  betas_cv1v2_pv2 <- sqrt(
                        Mod(
                          rbind(
                            do.call(rbind,lapply(Betas[group == 'C'               ],function(x) x[,pc])),
                            do.call(rbind,lapply(Betas[group == 'P' & visit == '2'],function(x) x[,pc]))
                          )
                        )
                      )
  labs_cv1v2_pv2 <- factor(c(rep('c_v1v2',sum(group=='C')),rep('p_v2',sum(group=='P' & visit=='2'))),levels=c('c_v1v2','p_v2'))
  ids_cv1v2_pv2 <- c(ident[group == 'C'],ident[group == 'P' & visit == '2'])

  betas_pv1_pv2 <- sqrt(
                      Mod(
                        rbind(
                          do.call(rbind,lapply(Betas[group == 'P' & visit == '1'],function(x) x[,pc])),
                          do.call(rbind,lapply(Betas[group == 'P' & visit == '2'],function(x) x[,pc]))
                        )
                      )
                    )
  labs_pv1_pv2 <- factor(c(rep('p_v1',sum(group=='P' & visit=='1')),rep('p_v2',sum(group=='P' & visit=='2'))),levels=c('p_v1','p_v2'))
#   ids_pv1_pv2 <- c(ident[group == 'P' & visit == '1'],ident[group == 'P' & visit == '2']) # not necessary

#   # student test p-values
#   pval_cv1v2_pv1 <- apply(betas_cv1v2_pv1,2,wrap_t.test,v=labs_cv1v2_pv1)
#   pval_cv1v2_pv2 <- apply(betas_cv1v2_pv2,2,wrap_t.test,v=labs_cv1v2_pv2)
#   pval_pv1_pv2 <- apply(betas_pv1_pv2,2,wrap_paired_t.test,v=labs_pv1_pv2)

  # wilcoxon test p-values
  pval_cv1v2_pv1 <- apply(betas_cv1v2_pv1,2,wrap_wilcox.test,v=labs_cv1v2_pv1)
  pval_cv1v2_pv2 <- apply(betas_cv1v2_pv2,2,wrap_wilcox.test,v=labs_cv1v2_pv2)
  pval_pv1_pv2 <- apply(betas_pv1_pv2,2,wrap_paired_wilcox.test,v=labs_pv1_pv2)

  plan(multicore,workers=32)

  # permutation test p-values as requested by JH
  nsim <- 1e3-1

  avgdiff_cv1v2_pv1 <- apply(betas_cv1v2_pv1,2,function(x,y) diff(tapply(x,y,mean)),y=labs_cv1v2_pv1)
  boot_avgdiff_cv1v2_pv1 <- listenv()
  for(sim in 1:nsim){
    boot_avgdiff_cv1v2_pv1[[sim]] %<-% {apply(betas_cv1v2_pv1,2,function(x,y) diff(tapply(x,y,mean)),y=permute_rep_meas(labs_cv1v2_pv1,ids_cv1v2_pv1))} %seed% T
  }
  boot_avgdiff_cv1v2_pv1 <- do.call(rbind,as.list(boot_avgdiff_cv1v2_pv1))
  pval_perm_cv1v2_pv1 <- rep(NA,length(avgdiff_cv1v2_pv1))
  for(col in 1:ncol(boot_avgdiff_cv1v2_pv1)){
    pval_perm_cv1v2_pv1[col] <- (1+sum(abs(boot_avgdiff_cv1v2_pv1[,col]) >= abs(avgdiff_cv1v2_pv1[col])))/(nsim+1)
  }

  avgdiff_cv1v2_pv2 <- apply(betas_cv1v2_pv2,2,function(x,y) diff(tapply(x,y,mean)),y=labs_cv1v2_pv2)
  boot_avgdiff_cv1v2_pv2 <- listenv()
  for(sim in 1:nsim){
    boot_avgdiff_cv1v2_pv2[[sim]] %<-% {apply(betas_cv1v2_pv2,2,function(x,y) diff(tapply(x,y,mean)),y=permute_rep_meas(labs_cv1v2_pv2,ids_cv1v2_pv2))} %seed% T
  }
  boot_avgdiff_cv1v2_pv2 <- do.call(rbind,as.list(boot_avgdiff_cv1v2_pv2))
  pval_perm_cv1v2_pv2 <- rep(NA,length(avgdiff_cv1v2_pv2))
  for(col in 1:ncol(boot_avgdiff_cv1v2_pv2)){
    pval_perm_cv1v2_pv2[col] <- (1+sum(abs(boot_avgdiff_cv1v2_pv2[,col]) >= abs(avgdiff_cv1v2_pv2[col])))/(nsim+1)
  }

  avgdiff_pv1_pv2 <- apply(betas_pv1_pv2,2,function(x,y) diff(tapply(x,y,mean)),y=labs_pv1_pv2)
  boot_avgdiff_pv1_pv2 <- listenv()
  for(sim in 1:nsim){
    boot_avgdiff_pv1_pv2[[sim]] %<-% {apply(betas_pv1_pv2,2,function(x,y) diff(tapply(x,y,mean)),y=sample(labs_pv1_pv2))} %seed% T
  }
  boot_avgdiff_pv1_pv2 <- do.call(rbind,as.list(boot_avgdiff_pv1_pv2))
  pval_perm_pv1_pv2 <- rep(NA,length(avgdiff_pv1_pv2))
  for(col in 1:ncol(boot_avgdiff_pv1_pv2)){
    pval_perm_pv1_pv2[col] <- (1+sum(abs(boot_avgdiff_pv1_pv2[,col]) >= abs(avgdiff_pv1_pv2[col])))/(nsim+1)
  }

  # histograms of p-values
  png(filename=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_hist_pval_cpc_',pc,'.png'),
      width=12,height=9,units='in',res=300)

  par(mfcol=c(3,2))

  hist(x=pval_cv1v2_pv1,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': HC x P V1'),sub='Student\'s t-test',xlab='p-value',col='gray66')

  hist(x=pval_cv1v2_pv2,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': HC x P V2'),sub='Student\'s t-test',xlab='p-value',col='gray66')

  hist(x=pval_pv1_pv2,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': P V1 x P V2'),sub='Student\'s t-test',xlab='p-value',col='gray66')

  hist(x=pval_perm_cv1v2_pv1,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': HC x P V1'),sub='Permutation test',xlab='p-value',col='gray66')

  hist(x=pval_perm_cv1v2_pv2,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': HC x P V2'),sub='Permutation test',xlab='p-value',col='gray66')

  hist(x=pval_perm_pv1_pv2,
       breaks=seq(from=0,to=1,by=0.01),
       main=paste0('cPC ',pc,': P V1 x P V2'),sub='Permutation test',xlab='p-value',col='gray66')
  dev.off()




  # cluster size permutation test
  nsim <- 1e3
  minsize <- 5
  alpha <- 0.005

  clustleft_cv1v2_pv1  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(pval_cv1v2_pv1[ind_left ] < alpha)],
                                               neigh=neigh_left,
                                               minsize=minsize,
                                               verbose=F)
  clustright_cv1v2_pv1 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(pval_cv1v2_pv1[ind_right] < alpha)],
                                               neigh=neigh_right,
                                               minsize=minsize,
                                               verbose=F)
  clustleft_cv1v2_pv2  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(pval_cv1v2_pv2[ind_left ] < alpha)],
                                               neigh=neigh_left,
                                               minsize=minsize,
                                               verbose=F)
  clustright_cv1v2_pv2 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(pval_cv1v2_pv2[ind_right] < alpha)],
                                               neigh=neigh_right,
                                               minsize=minsize,
                                               verbose=F)
  clustleft_pv1_pv2  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(pval_pv1_pv2[ind_left ] < alpha)],
                                             neigh=neigh_left,
                                             minsize=minsize,
                                             verbose=F)
  clustright_pv1_pv2 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(pval_pv1_pv2[ind_right] < alpha)],
                                             neigh=neigh_right,
                                             minsize=minsize,
                                             verbose=F)

  # HC x P V1

  boot_maxclustsize_cv1v2_pv1 <- listenv()

  for(sim in 1:nsim){

    boot_maxclustsize_cv1v2_pv1[[sim]] %<-% {
#       boot_pval_cv1v2_pv1 <- apply(betas_cv1v2_pv1,2,wrap_t.test,v=permute_rep_meas(labs_cv1v2_pv1,ids_cv1v2_pv1)) # student t-test
      boot_pval_cv1v2_pv1 <- apply(betas_cv1v2_pv1,2,wrap_wilcox.test,v=permute_rep_meas(labs_cv1v2_pv1,ids_cv1v2_pv1)) # wilcoxon test
      boot_clustleft_cv1v2_pv1  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(boot_pval_cv1v2_pv1[ind_left ] < alpha)],
                                                        neigh=neigh_left,
                                                        verbose=F)
      boot_clustright_cv1v2_pv1 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(boot_pval_cv1v2_pv1[ind_right] < alpha)],
                                                        neigh=neigh_right,
                                                        verbose=F)

      boot_maxclustleft_cv1v2_pv1 <- NA
      if(length(boot_clustleft_cv1v2_pv1) > 0){boot_maxclustleft_cv1v2_pv1 <- max(sapply(boot_clustleft_cv1v2_pv1,function(x) x$size))}
      boot_maxclustright_cv1v2_pv1 <- NA
      if(length(boot_clustright_cv1v2_pv1) > 0){boot_maxclustright_cv1v2_pv1 <- max(sapply(boot_clustright_cv1v2_pv1,function(x) x$size))}

      if(is.na(boot_maxclustleft_cv1v2_pv1) & is.na(boot_maxclustright_cv1v2_pv1)){
        return(0)
      }else{
        return(max(c(boot_maxclustleft_cv1v2_pv1,boot_maxclustright_cv1v2_pv1),na.rm=T))
      }

    } %seed% T
  }

  boot_maxclustsize_cv1v2_pv1 <- do.call(c,as.list(boot_maxclustsize_cv1v2_pv1))


  cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))
  cat(paste0('cPC',pc,': HC x P V1\n'))
  cat(paste0('Number of clusters in left hemisphere: ',length(clustleft_cv1v2_pv1),'\n'))
  cat('Cluster sizes in left hemisphere:\n')
  cat(sort(sapply(clustleft_cv1v2_pv1,function(x) x$size),decreasing=T),'\n')
  cat(paste0('Number of clusters in right hemisphere: ',length(clustright_cv1v2_pv1),'\n'))
  cat('Cluster sizes in right hemisphere:\n')
  cat(sort(sapply(clustright_cv1v2_pv1,function(x) x$size),decreasing=T),'\n')
  cat(paste0('95% quantile of the maximum cluster size distribution: ',quantile(boot_maxclustsize_cv1v2_pv1,probs=0.95,na.rm=T),'\n'))
  cat(paste0('Number of clusters in the left hemisphere with size above the limit: ',sum(sapply(clustleft_cv1v2_pv1,function(x) x$size) > quantile(boot_maxclustsize_cv1v2_pv1,probs=0.95,na.rm=T)),'\n'))
  cat(paste0('Number of clusters in the right hemisphere with size above the limit: ',sum(sapply(clustright_cv1v2_pv1,function(x) x$size) > quantile(boot_maxclustsize_cv1v2_pv1,probs=0.95,na.rm=T)),'\n'))
  cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))


  boot_maxclustsize_cv1v2_pv2 <- listenv()

  for(sim in 1:nsim){

    boot_maxclustsize_cv1v2_pv2[[sim]] %<-% {

#       boot_pval_cv1v2_pv2 <- apply(betas_cv1v2_pv2,2,wrap_t.test,v=permute_rep_meas(labs_cv1v2_pv2,ids_cv1v2_pv2)) # students t-test
      boot_pval_cv1v2_pv2 <- apply(betas_cv1v2_pv2,2,wrap_wilcox.test,v=permute_rep_meas(labs_cv1v2_pv2,ids_cv1v2_pv2)) # wilcoxon test
      boot_clustleft_cv1v2_pv2  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(boot_pval_cv1v2_pv2[ind_left ] < alpha)],
                                                        neigh=neigh_left,
                                                        verbose=F)
      boot_clustright_cv1v2_pv2 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(boot_pval_cv1v2_pv2[ind_right] < alpha)],
                                                        neigh=neigh_right,
                                                        verbose=F)

      boot_maxclustleft_cv1v2_pv2 <- NA
      if(length(boot_clustleft_cv1v2_pv2) > 0){boot_maxclustleft_cv1v2_pv2 <- max(sapply(boot_clustleft_cv1v2_pv2,function(x) x$size))}
      boot_maxclustright_cv1v2_pv2 <- NA
      if(length(boot_clustright_cv1v2_pv2) > 0){boot_maxclustright_cv1v2_pv2 <- max(sapply(boot_clustright_cv1v2_pv2,function(x) x$size))}

      if(is.na(boot_maxclustleft_cv1v2_pv2) & is.na(boot_maxclustright_cv1v2_pv2)){
        return(0)
      }else{
        return(max(c(boot_maxclustleft_cv1v2_pv2,boot_maxclustright_cv1v2_pv2),na.rm=T))
      }

    } %seed% T
  }

  boot_maxclustsize_cv1v2_pv2 <- do.call(c,as.list(boot_maxclustsize_cv1v2_pv2))


  cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))
  cat(paste0('cPC',pc,': HC x P V2\n'))
  cat(paste0('Number of clusters in left hemisphere: ',length(clustleft_cv1v2_pv2),'\n'))
  cat('Cluster sizes in left hemisphere:\n')
  cat(sort(sapply(clustleft_cv1v2_pv2,function(x) x$size),decreasing=T),'\n')
  cat(paste0('Number of clusters in right hemisphere: ',length(clustright_cv1v2_pv2),'\n'))
  cat('Cluster sizes in right hemisphere:\n')
  cat(sort(sapply(clustright_cv1v2_pv2,function(x) x$size),decreasing=T),'\n')
  cat(paste0('95% quantile of the maximum cluster size distribution: ',quantile(boot_maxclustsize_cv1v2_pv2,probs=0.95,na.rm=T),'\n'))
  cat(paste0('Number of clusters in the left hemisphere with size above the limit: ',sum(sapply(clustleft_cv1v2_pv2,function(x) x$size) > quantile(boot_maxclustsize_cv1v2_pv2,probs=0.95,na.rm=T)),'\n'))
  cat(paste0('Number of clusters in the right hemisphere with size above the limit: ',sum(sapply(clustright_cv1v2_pv2,function(x) x$size) > quantile(boot_maxclustsize_cv1v2_pv2,probs=0.95,na.rm=T)),'\n'))
  cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))

#   boot_maxclustsize_pv1_pv2 <- listenv()
#
#   for(sim in 1:nsim){
#
#     boot_maxclustsize_pv1_pv2[[sim]] %<-% {
#
# #       boot_pval_pv1_pv2 <- apply(betas_pv1_pv2,2,wrap_paired_t.test,v=sample(labs_pv1_pv2)) # students t-test
#       boot_pval_pv1_pv2 <- apply(betas_pv1_pv2,2,wrap_paired_wilcox.test,v=sample(labs_pv1_pv2)) # wilcoxon test
#       boot_clustleft_pv1_pv2  <- compute_cluster_size(sigind=(hemi_left[veri_left])[which(boot_pval_pv1_pv2[ind_left ] < alpha)],
#                                                       neigh=neigh_left,
#                                                       verbose=F)
#       boot_clustright_pv1_pv2 <- compute_cluster_size(sigind=(hemi_right[veri_right] - tail(hemi_left,1))[which(boot_pval_pv1_pv2[ind_right] < alpha)],
#                                                       neigh=neigh_right,
#                                                       verbose=F)
#
#       boot_maxclustleft_pv1_pv2 <- NA
#       if(length(boot_clustleft_pv1_pv2) > 0){boot_maxclustleft_pv1_pv2 <- max(sapply(boot_clustleft_pv1_pv2,function(x) x$size))}
#       boot_maxclustright_pv1_pv2 <- NA
#       if(length(boot_clustright_pv1_pv2) > 0){boot_maxclustright_pv1_pv2 <- max(sapply(boot_clustright_pv1_pv2,function(x) x$size))}
#
#       if(is.na(boot_maxclustleft_pv1_pv2) & is.na(boot_maxclustright_pv1_pv2)){
#         return(0)
#       }else{
#         return(max(c(boot_maxclustleft_pv1_pv2,boot_maxclustright_pv1_pv2),na.rm=T))
#       }
#
#     } %seed% T
#   }
#
#   boot_maxclustsize_pv1_pv2 <- do.call(c,as.list(boot_maxclustsize_pv1_pv2))
#
#
#   cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))
#   cat(paste0('cPC',pc,': P V1 x P V2\n'))
#   cat(paste0('Number of clusters in left hemisphere: ',length(clustleft_pv1_pv2),'\n'))
#   cat('Cluster sizes in left hemisphere:\n')
#   cat(sort(sapply(clustleft_pv1_pv2,function(x) x$size),decreasing=T),'\n')
#   cat(paste0('Number of clusters in right hemisphere: ',length(clustright_pv1_pv2),'\n'))
#   cat('Cluster sizes in right hemisphere:\n')
#   cat(sort(sapply(clustright_pv1_pv2,function(x) x$size),decreasing=T),'\n')
#   cat(paste0('95% quantile of the maximum cluster size distribution: ',quantile(boot_maxclustsize_pv1_pv2,probs=0.95,na.rm=T),'\n'))
#   cat(paste0('Number of clusters in the left hemisphere with size above the limit: ',sum(sapply(clustleft_pv1_pv2,function(x) x$size) > quantile(boot_maxclustsize_pv1_pv2,probs=0.95,na.rm=T)),'\n'))
#   cat(paste0('Number of clusters in the right hemisphere with size above the limit: ',sum(sapply(clustright_pv1_pv2,function(x) x$size) > quantile(boot_maxclustsize_pv1_pv2,probs=0.95,na.rm=T)),'\n'))
#   cat(paste0(paste0(rep('#',32),collapse=' '),'\n'))

  Res[[pc]] <- list(
                 hcv1v2_pv1=list(nsim=nsim,
                                 minsize=minsize,
                                 clust_left=clustleft_cv1v2_pv1,
                                 clust_right=clustright_cv1v2_pv1,
                                 boot_maxclustsize=boot_maxclustsize_cv1v2_pv1),
                 hcv1v2_pv2=list(nsim=nsim,
                                 minsize=minsize,
                                 clust_left=clustleft_cv1v2_pv2,
                                 clust_right=clustright_cv1v2_pv2,
                                 boot_maxclustsize=boot_maxclustsize_cv1v2_pv2)) #,
#                  pv1_pv2=list(nsim=nsim,
#                               minsize=minsize,
#                               clust_left=clustleft_pv1_pv2,
#                               clust_right=clustright_pv1_pv2,
#                               boot_maxclustsize=boot_maxclustsize_pv1_pv2))

  # simple averages for visualizations
  avg_hc_v1v2 <- apply(do.call(rbind,lapply(Betas[group == 'C'               ],function(x) x[,pc])),2,mean)
  avg_p__v1   <- apply(do.call(rbind,lapply(Betas[group == 'P' & visit == '1'],function(x) x[,pc])),2,mean)
  avg_p__v2   <- apply(do.call(rbind,lapply(Betas[group == 'P' & visit == '2'],function(x) x[,pc])),2,mean)

  avgdiff_hc_v1v2_p_v1 <- avg_p__v1 - avg_hc_v1v2
  avgdiff_hc_v1v2_p_v2 <- avg_p__v2 - avg_hc_v1v2

  # visualizing the modulus
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- sqrt(Mod(avgdiff_hc_v1v2_p_v1[ind_left]))
  visualize_on_surface(x=helper,
                      surf=inflated_left,
                      file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_cpc',pc,'_hc_v1v2_p_v1.html'))

  # visualizing the modulus
  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- sqrt(Mod(avgdiff_hc_v1v2_p_v1[ind_right]))
  visualize_on_surface(x=helper,
                      surf=inflated_right,
                      file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_cpc',pc,'_hc_v1v2_p_v1.html'))

  # visualizing the modulus
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- sqrt(Mod(avgdiff_hc_v1v2_p_v2[ind_left]))
  visualize_on_surface(x=helper,
                      surf=inflated_left,
                      file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_cpc',pc,'_hc_v1v2_p_v2.html'))

  # visualizing the modulus
  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- sqrt(Mod(avgdiff_hc_v1v2_p_v2[ind_right]))
  visualize_on_surface(x=helper,
                      surf=inflated_right,
                      file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_cpc',pc,'_hc_v1v2_p_v2.html'))

  # P - VALUES
  # C V1v2 x P V1
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_cv1v2_pv1[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_cpc',pc,'_hc_v1v2_p_v1.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_cv1v2_pv1[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_cpc',pc,'_hc_v1v2_p_v1.html'))

  # C V1V2 x P V2
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_cv1v2_pv2[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_cpc',pc,'_hc_v1v2_p_v2.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_cv1v2_pv2[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_cpc',pc,'_hc_v1v2_p_v2.html'))

  # P V1 x P V2
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_pv1_pv2[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_cpc',pc,'_p_v1_p_v2.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_pv1_pv2[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_cpc',pc,'_p_v1_p_v2.html'))



  # PERMUTATION TEST P - VALUES
  # C V1v2 x P V1
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_perm_cv1v2_pv1[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_perm_cpc',pc,'_hc_v1v2_p_v1.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_perm_cv1v2_pv1[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_perm_cpc',pc,'_hc_v1v2_p_v1.html'))

  # C V1V2 x P V2
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_perm_cv1v2_pv2[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_perm_cpc',pc,'_hc_v1v2_p_v2.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_perm_cv1v2_pv2[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_perm_cpc',pc,'_hc_v1v2_p_v2.html'))

  # P V1 x P V2
  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- pval_perm_pv1_pv2[ind_left]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_left,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_left_mod_pval_perm_cpc',pc,'_p_v1_p_v2.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- pval_perm_pv1_pv2[ind_right]
  cols <- pval_sig_cols(helper)

  visualize_on_surface2(x=helper,
                        surf=inflated_right,
                        cols=cols,
                        file=paste0(path_rend,'/',foldname_rend_comp,'/',dec,'_',smooth,'_cpca_global_pimvr_right_mod_pval_perm_cpc',pc,'_p_v1_p_v2.html'))
}

plan(sequential)

# save('Res',file=paste0(path_work,'workspace_',dec,'_',smooth,'_cluster_size_inference_',gsub(x=Sys.Date(),pattern='-',replacement=''),'.Rdata'))
