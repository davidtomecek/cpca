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

dec <- 'fsaverage5'
fs <- 1

path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')
path_labe <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/label/')

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'
# file_work <- 'workspace_fsaverage5_cpca_pimvr_20251001.RData'
file_work <- 'workspace_fsaverage5_cpca_pimvr_partial_jd.RData'

path_rend <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/'
foldname_rend <- 'cpa_global_pimvr_fsaverage5_jd'
if(!dir.exists(paste0(path_rend,'/',foldname_rend,'/'))){dir.create(paste0(path_rend,'/',foldname_rend,'/'),recursive=T)}

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
# ncl <- 128
# cmap <- hcl.colors(n=ncl,palette='Zissou 1')

# loading the precomputed group cPCA data
# load(paste0(path_work,'/fsaverage6_cpca_pimvr_check_20251104.RData'))
# load(paste0(path_work,'/fsaverage6_cpca_pimvr_check_help_20251104.RData'))

# load(paste0(path_work,'workspace_fsaverage6_sm5_cpca_pimvr_partial.RData'))
#
# load(paste0(path_work,'ESO_P03681_20250602_0850_2/fsaverage6_sm5_group_cpca_preprocessed_data.RData'))

load(paste0(path_work,file_work))
load(paste0(path_work,'ESO_P03681_20250602_0850_2/fsaverage5_group_cpca_preprocessed_data.RData'))

size_help <- nrow(data_preproc)

n_pc <- ncol(cPCA_list[[1]]$W)

#     # select the result with the greatest overall explained variance
#     R2 <- do.call(rbind,lapply(cPCA_list,function(x) x$R2))
#     R2_ord <- apply(R2,2,order,decreasing=T)
#
#     ind_best <- table(R2_ord[1,])
#     ind_best <- as.integer(names(ind_best)[which.max(ind_best)])

    W_glob <- cPCA_list[[1]]$W

    W_glob_adj <- matrix(NA,nrow=size_help,ncol=n_pc)
    W_glob_adj[-coni,] <- W_glob

# # preparing ranges for the different cPCs
# pc_ranges <- data.frame(pc=1:n_pc,ll_mod=rep(NA,times=n_pc),ul_mod=rep(NA,times=n_pc))
# for(pc in 1:n_pc){
#   ll <- sapply(CPCA,function(x) min(sqrt(Mod(x$x[c(ind_left,ind_right),pc])),na.rm=T))
#   ul <- sapply(CPCA,function(x) max(sqrt(Mod(x$x[c(ind_left,ind_right),pc])),na.rm=T))
#
#   pc_ranges$ll_mod[pc] <- min(ll,na.rm=T)
#   pc_ranges$ul_mod[pc] <- max(ul,na.rm=T)
# }

# going through the subsets
for(set in 1:length(Sets)){

#   # scree plot
#   png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca4set_opt_',names(Sets)[i],'_expvar.png'),
#       width=4,height=4,units='in',res=300)
#   plot(CPCA[[i]]$sdev^2/sum(CPCA[[i]]$sdev^2),xlim=c(1,10),main=paste0('Group ',names(Sets)[i]),xlab='Principal component',ylab='Explained variance [%]',type='b',frame=F)
#   dev.off()

#   # loadings plot - this should be divided by datasets and compared to the individual results
#   for(fold in 1:length(Sets[[set]])){
#
#
#     load(paste0(path_work,Sets[[set]][fold],'/',dec,'_group_cpca_preprocessed_data.RData'))
#     data_preproc <- t(data_preproc[-coni,])
#     scores_indi <- data_preproc %*% W_glob
#
#     png(filename=paste0(path_rend,foldname_rend,'/',dec,'_cpca_global_pimvr_data_',Sets[[set]][fold],'_scores.png'),
#         width=8,height=2*n_pc,units='in',res=300)
#     par(mfrow=c(n_pc,2))
#     for(pc in 1:n_pc){
#       plot(x=seq(from=0,by=1/fs,length.out=Sizes[[fold]][2]),
#            y=Mod(scores_indi[,pc]),
#            xlab='Time [s]',ylab=paste0('Scores ',pc,' - Modulus'),type='l',frame=F,lwd=1.2)
#
#         plot(x=seq(from=0,by=1/fs,length.out=Sizes[[fold]][2]),
#              y=Arg(scores_indi[,pc]),
#              xlab='Time [s]',ylab=paste0('Scores ',pc,' - Argument'),type='l',frame=F,lwd=1.2)
#       }
#       dev.off()
#   }

  for(pc in 1:n_pc){

    phase_adj <- adjust_phase(Arg(W_glob_adj[,pc]))

    # left hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- phase_adj[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca_global_pimvr_left_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- (Mod(W_glob_adj[ind_left,pc]))
    visualize_on_surface(x=helper,
                         surf=inflated_left,
#                          range=range(Mod(W_glob[,pc])), # *c(1.1,0.9), # c(pc_ranges$ll_mod[pc],pc_ranges$ul_mod[pc]),
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca_global_pimvr_left_mod_cpc_',pc,'.html'))

    # right hemisphere
    # visualizing the phase
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- phase_adj[ind_right]
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca_global_pimvr_right_arg_cpc_',pc,'.html'))

    # visualizing the modulus
    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- (Mod(W_glob_adj[ind_right,pc]))
    visualize_on_surface(x=helper,
                         surf=inflated_right,
#                          range=range(Mod(W_glob[,pc])), #*c(1.1,0.9), # c(pc_ranges$ll_mod[pc],pc_ranges$ul_mod[pc]),
                         file=paste0(path_rend,'/',foldname_rend,'/',dec,'_cpca_global_pimvr_right_mod_cpc_',pc,'.html'))

  }

}
