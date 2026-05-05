# TODO
# perhaps a summary visualization of all the experiments - variance in scores across all experiments
# some conclusions



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
    ind <- round((ncols-1)*(x - range[1])/(range[2] - range[1])+1)
    ind[ind <= 1] <- 1
    ind[ind > ncols] <- ncols
    cols <- cmap[ind]
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

  unlink(x=paste0(strsplit(x=file,split='.',fixed=T)[[1]][1],'_files'))

}

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

files_cpca_orig <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9]{8}\\.RData$'))
files_cpca_expe <- list.files(path_work,pattern=paste0('^workspace_cpca4set_magnexpe_',dec,'_[0-9]{8}\\.RData$'))

set <- 1

load(paste0(path_work,tail(files_cpca_orig,n=1)))
rm('Data_comb')
for(deleset in (1:length(CPCA))[-set]){CPCA[[deleset]] <- NULL}
gc()

load(paste0(path_work,tail(files_cpca_expe,n=1)))

path_imag <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20241122/cpca4set/magnexpe/'
if(!dir.exists(path_imag)){dir.create(path_imag,recursive=T)}

path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')

library(freesurferformats)
inflated_left  <- read.fs.surface(paste0(path_surf,'lh.inflated'))
inflated_right <- read.fs.surface(paste0(path_surf,'rh.inflated'))

ind_left  <- 1:(length(veri_left))
ind_right <- (1+length(veri_left)):(length(veri_left) + length(veri_right))

expe_par <- data.frame(n=sapply(CPCA_magnexpe[[set]],function(x) x$expe_n),
                       ind=sapply(CPCA_magnexpe[[set]],function(x) paste0('(',paste0(x$expe_ind,collapse=', '),')')),
                       coef=sapply(CPCA_magnexpe[[set]],function(x) paste0('(',paste0(x$expe_coef,collapse=', '),')')))

expe_ntry <- length(CPCA_magnexpe[[1]])

cmap <- rainbow(n=expe_ntry)

nr <- ceiling(sqrt(expe_ntry))
nc <- ceiling(expe_ntry/nr)

fs <- 1
n_dataset <- length(Sets[[set]])
n_samp <- nrow(CPCA[[set]]$rotation)/n_dataset
n_pc <- 3

for(pc in 1:n_pc){

  # visualizaiton of loadings

  # one long time-series

  expe_range <- range(rbind(range(Mod(CPCA[[set]]$rotation[,pc])),
                            do.call(rbind,lapply(CPCA_magnexpe[[set]],function(x) range(Mod(x$expe_cpca$rotation[,pc]))))))


  png(filename=paste0(path_imag,'loadings_long_pc',formatC(pc,digits=2,flag='0'),'.png'),
      width=n_dataset*4,height=3,res=300,units='in')

  plot(x=seq(from=0,by=1/fs,length.out=nrow(CPCA[[set]]$rotation)),
       xlab='Time [s]',
       y=Mod(CPCA[[set]]$rotation[,pc]),
       ylab=paste0('Mod of cPC',pc,' loading []'),
       ylim=expe_range,
       main='',
       type='l',lwd=1.2,col='black',frame=F)
    for(expe in 1:expe_ntry){
      lines(x=seq(from=0,by=1/fs,length.out=nrow(CPCA[[set]]$rotation)),
            y=Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$rotation[,pc]),
            col=cmap[expe])
    }
  dev.off()

  # divided per dataset

  png(filename=paste0(path_imag,'loadings_divided_pc',formatC(pc,digits=2,flag='0'),'.png'),
      width=nc*4,height=nr*3,res=300,units='in')

  par(mfrow=c(nr,nc))

  for(dataset in 1:length(Sets[[set]])){

    plot(x=seq(from=0,by=1/fs,length.out=n_samp),
         xlab='Time [s]',
         y=Mod(CPCA[[set]]$rotation[(n_samp*(dataset-1)+1):(n_samp*dataset),pc]),
         ylab=paste0('Mod of cPC',pc,' loading []'),
         ylim=expe_range,
         main=Sets[[set]][dataset],
         type='l',lwd=1.2,col='black',frame=F)
    for(expe in 1:expe_ntry){
      lines(x=seq(from=0,by=1/fs,length.out=n_samp),
            y=Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$rotation[(n_samp*(dataset-1)+1):(n_samp*dataset),pc]),
            col=cmap[expe])
    }
  }
  dev.off()

  # visualization of scores

  for(expe in 1:expe_ntry){
    # modulus
    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$x[ind_left,pc])
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         file=paste0(path_imag,'scores_mod_left_expe',expe,'_cpc',pc,'.html'))

    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$x[ind_right,pc])
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         file=paste0(path_imag,'scores_mod_right_expe',expe,'_cpc',pc,'.html'))


    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- Mod(CPCA[[set]]$x[ind_left,pc]) - Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$x[ind_left,pc])
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         file=paste0(path_imag,'scores_mod_diff_left_expe',expe,'_cpc',pc,'.html'),
                         range=c(-1,1)*max(Mod(CPCA[[set]]$x[,pc])))

    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- Mod(CPCA[[set]]$x[ind_right,pc]) - Mod(CPCA_magnexpe[[set]][[expe]]$expe_cpca$x[ind_right,pc])
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         file=paste0(path_imag,'scores_mod_diff_right_expe',expe,'_cpc',pc,'.html'),
                         range=c(-1,1)*max(Mod(CPCA[[set]]$x[,pc])))

    # phase
    phase_orig_adj <- adjust_phase(Arg(CPCA[[set]]$x[,pc]))
    phase_expe_adj <- adjust_phase(Arg(CPCA_magnexpe[[set]][[expe]]$expe_cpca$x[,pc]))

    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- phase_expe_adj[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         file=paste0(path_imag,'scores_arg_left_expe',expe,'_cpc',pc,'.html'))

    helper <- rep(NA,times=nrow(inflated_right$vertices))
    helper[hemi_right[veri_right] - max(hemi_left)] <- phase_expe_adj[ind_right]
    visualize_on_surface(x=helper,
                         surf=inflated_right,
                         file=paste0(path_imag,'scores_arg_right_expe',expe,'_cpc',pc,'.html'))


    helper <- rep(NA,times=nrow(inflated_left$vertices))
    helper[hemi_left[veri_left]] <- phase_orig_adj[ind_left] - phase_expe_adj[ind_left]
    visualize_on_surface(x=helper,
                         surf=inflated_left,
                         file=paste0(path_imag,'scores_arg_diff_left_expe',expe,'_cpc',pc,'.html'))

     helper <- rep(NA,times=nrow(inflated_right$vertices))
     helper[hemi_right[veri_right] - max(hemi_left)] <- phase_orig_adj[ind_right] - phase_expe_adj[ind_right]
     visualize_on_surface(x=helper,
                          surf=inflated_right,
                          file=paste0(path_imag,'scores_arg_diff_right_expe',expe,'_cpc',pc,'.html'))
  }

#   for(pc in 1:n_pc){
  # summary of all experiments - variance in scores

  library(complexlm)

  expe_sumvar <- apply(do.call(cbind,lapply(CPCA_magnexpe[[set]],function(x) x$expe_cpca$x[,pc])),1,var)

  helper <- rep(NA,times=nrow(inflated_left$vertices))
  helper[hemi_left[veri_left]] <- expe_sumvar[ind_left]
  visualize_on_surface(x=helper,
                       surf=inflated_left,
                       file=paste0(path_imag,'scores_sumvar_left_cpc',pc,'.html'))

  helper <- rep(NA,times=nrow(inflated_right$vertices))
  helper[hemi_right[veri_right] - max(hemi_left)] <- expe_sumvar[ind_right]
  visualize_on_surface(x=helper,
                       surf=inflated_right,
                       file=paste0(path_imag,'scores_sumvar_right_cpc',pc,'.html'))

}
