rm(list=ls(all.names=T))
graphics.off()

Sys.setenv(RGL_USE_NULL=T)

adjust_phase_cost <- function(p,x){
  x <- ((x - p) %% (2*pi)) - pi
  abs(mean(x))
}

adjust_phase <- function(p){
  o <- optim(par=0,fn=adjust_phase_cost,x=p,method='Brent',lower=-pi,upper=pi)$par
  q <- ((p - o) %% (2*pi)) - pi
  return(q)
}

visualize_on_surface <- function(x,surf,file,cmap='Zissou 1',ncols=256,range=NULL){

  cmap <- hcl.colors(n=ncols,palette=cmap)

  if(is.null(range)){
    cols <- cmap[round((ncols-1)*(x - min(x,na.rm=T))/(max(x,na.rm=T) - min(x,na.rm=T))+1)]
    cols[is.na(x)] <- 'gray66'
  }else{
    cols <- cmap[round((ncols-1)*(x - range[1])/(range[2] - range[1])+1)]
    cols[is.na(x)] <- 'gray66'
  }

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

dec <- 'fsaverage6'
n_pc <- 3
fs_group <- 2
fs_indiv <- 3

path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
# file_work <- list.files(path_work,pattern=paste0('^workspace_cpca4set_',dec,'_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\\.RData'))
# file_work <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\\.RData'))
file_work <- list.files(path_work,pattern=paste0('^workspace_cpca4set_optnorm_cpca_',dec,'_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]\\.RData'))
file_work <- tail(file_work,n=1)

path_rend <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20241122/'
foldname_rend <- 'cpca4set_optnorm'
if(!dir.exists(paste0(path_rend,'/',foldname_rend,'/'))){dir.create(paste0(path_rend,'/',foldname_rend,'/'),recursive=T)}

library(freesurferformats)
inflated_left  <- read.fs.surface(paste0(path_surf,'lh.inflated'))
inflated_right <- read.fs.surface(paste0(path_surf,'rh.inflated'))

# ncl <- 128
# cmap <- hcl.colors(n=ncl,palette='Zissou 1')

# loading the precomputed group cPCA data
load(paste0(path_work,'/',file_work))

# removing some loadings and scores to reduce the memory load
for(i in 1:length(CPCA)){
  CPCA[[i]]$rotation <- CPCA[[i]]$rotation[,1:n_pc]
  CPCA[[i]]$center <- NULL
  CPCA[[i]]$scale <- NULL
  CPCA[[i]]$x <- CPCA[[i]]$x[,1:n_pc]
}
gc(); gc(); gc()

ind_left  <- 1:(length(veri_left))
ind_right <- (1+length(veri_left)):(length(veri_left) + length(veri_right))

# preparing ranges for the different cPCs
pc_ranges <- data.frame(pc=1:n_pc,ll_mod=rep(NA,times=n_pc),ul_mod=rep(NA,times=n_pc))
for(pc in 1:n_pc){
  ll <- sapply(CPCA,function(x) min(sqrt(Mod(x$x[c(ind_left,ind_right),pc])),na.rm=T))
  ul <- sapply(CPCA,function(x) max(sqrt(Mod(x$x[c(ind_left,ind_right),pc])),na.rm=T))

  pc_ranges$ll_mod[pc] <- min(ll,na.rm=T)
  pc_ranges$ul_mod[pc] <- max(ul,na.rm=T)
}

# going through the subsets
for(i in 1:length(CPCA)){

  # scree plot
#   png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_expvar.png'),
  png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_expvar.png'),
      width=4,height=4,units='in',res=300)
  plot(CPCA[[i]]$sdev^2/sum(CPCA[[i]]$sdev^2),xlim=c(1,10),main=paste0('Group ',names(Sets)[i]),xlab='Principal component',ylab='Explained variance [%]',type='b',frame=F)
  dev.off()

  # loadings plot - this should be divided by datasets and compared to the individual results
  for(j in 1:length(Sets[[i]])){

    try({

      load(paste0(path_work,Sets[[i]][j],'/',dec,'_cpca.RData'))
      cpca1$rotation <- cpca1$rotation[,1:n_pc]
      cpca1$center <- NULL
      cpca1$scale <- NULL
      cpca1$x <- NULL

      n_group <- length(CPCA[[i]]$rotation[,pc])/length(Sets[[i]])
      n_indiv <- length(cpca1$rotation[,pc])

#       png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_data_',Sets[[i]][j],'_loadings.png'),
#       png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_data_',Sets[[i]][j],'_loadings.png'),
      png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_optnorm_',names(Sets)[i],'_data_',Sets[[i]][j],'_loadings.png'),
          width=8,height=2*n_pc,units='in',res=300)
      par(mfrow=c(n_pc,2))
      for(pc in 1:n_pc){
        plot(x=seq(from=0,by=1/fs_group,length.out=n_group),y=Mod(CPCA[[i]]$rotation[((j-1)*n_group+1):(j*n_group),pc]),ylim=range(c(Mod(CPCA[[i]]$rotation[((j-1)*n_group+1):(j*n_group),pc]),Mod(cpca1$rotation[,pc]))),xlab='Time [s]',ylab=paste0('Loading ',pc,' - Modulus'),type='l',frame=F,lwd=1.2)
        lines(x=seq(from=0,by=1/fs_indiv,length.out=n_indiv),y=Mod(cpca1$rotation[,pc]),lwd=0.4)

        plot(x=seq(from=0,by=1/fs_group,length.out=n_group),y=Arg(CPCA[[i]]$rotation[((j-1)*n_group+1):(j*n_group),pc]),xlab='Time [s]',ylab=paste0('Loading ',pc,' - Argument'),type='l',frame=F,lwd=1.2)
        lines(x=seq(from=0,by=1/fs_indiv,length.out=n_indiv),y=Arg(cpca1$rotation[,pc]),lwd=0.4)
      }
      dev.off()
    },
    silent=T)
  }

  for(pc in 1:n_pc){

    phase_adj <- adjust_phase(Arg(CPCA[[i]]$x[,pc]))

    # left hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- phase_adj[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_left_arg_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_left_arg_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_',names(Sets)[i],'_left_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- sqrt(Mod(CPCA[[i]]$x[ind_left,pc]))
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         range=c(pc_ranges$ll_mod[pc],pc_ranges$ul_mod[pc]),
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_left_mod_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_left_mod_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_',names(Sets)[i],'_left_mod_cpc_',pc,'.html'))

    # right hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- phase_adj[ind_right]
    visualize_on_surface(x=helper,
                         surf=inflated_right,
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_right_arg_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_right_arg_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_',names(Sets)[i],'_right_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- sqrt(Mod(CPCA[[i]]$x[ind_right,pc]))
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         range=c(pc_ranges$ll_mod[pc],pc_ranges$ul_mod[pc]),
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_',names(Sets)[i],'_right_mod_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_right_mod_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_',names(Sets)[i],'_right_mod_cpc_',pc,'.html'))

  }

}

# TODO
# differences between the groups, perhaps gpv1 x gpv2, gcv1 x gcv2, gpv1 x gcv1 and gpv2 x gcv2

diff_list <- cbind(c('gpv1','gcv1','gpv1','gpv2'),
                   c('gpv2','gcv2','gcv1','gcv2'))

pc_ranges_diff <- data.frame(pc=1:n_pc,ll_mod=rep(NA,times=n_pc),ul=rep(NA,times=n_pc))

for(pc in 1:n_pc){

  ll <- rep(NA,times=nrow(diff_list))
  ul <- rep(NA,times=nrow(diff_list))

  for(i in 1:nrow(diff_list)){

    ind1 <- which(names(Sets) == diff_list[i,1])
    ind2 <- which(names(Sets) == diff_list[i,2])

    mod_diff <- Mod(CPCA[[ind2]]$x[,pc]) - Mod(CPCA[[ind1]]$x[,pc])
    mod_diff <- sign(mod_diff)*sqrt(abs(mod_diff))

    ll[i] <- min(mod_diff,na.rm=T)
    ul[i] <- max(mod_diff,na.rm=T)

  }

  pc_ranges_diff$ll_mod[pc] <- min(ll,na.rm=T)
  pc_ranges_diff$ul_mod[pc] <- max(ul,na.rm=T)

}

for(i in 1:nrow(diff_list)){

  ind1 <- which(names(Sets) == diff_list[i,1])
  ind2 <- which(names(Sets) == diff_list[i,2])

  for(pc in 1:n_pc){
    phase_diff_adj <- adjust_phase(Arg(CPCA[[ind2]]$x[,pc]) - Arg(CPCA[[ind1]]$x[,pc]))
    mod_diff <- Mod(CPCA[[ind2]]$x[,pc]) - Mod(CPCA[[ind1]]$x[,pc])
    mod_diff <- sign(mod_diff)*sqrt(abs(mod_diff))

    # left hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- phase_diff_adj[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_arg_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_arg_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- mod_diff[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         range=c(pc_ranges_diff$ll_mod[pc],pc_ranges_diff$ul_mod[pc]),
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_mod_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_mod_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_diff_',diff_list[i,2],'_',diff_list[i,1],'_left_mod_cpc_',pc,'.html'))

    # right hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- phase_diff_adj[ind_right]
    visualize_on_surface(x=helper,
                         surf=inflated_right,
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_arg_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_arg_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- mod_diff[ind_right]
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         range=c(pc_ranges_diff$ll_mod[pc],pc_ranges_diff$ul_mod[pc]),
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_mod_cpc_',pc,'.html'))
#                          file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_opt_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_mod_cpc_',pc,'.html'))
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca4set_optnorm_diff_',diff_list[i,2],'_',diff_list[i,1],'_right_mod_cpc_',pc,'.html'))
  }
}
