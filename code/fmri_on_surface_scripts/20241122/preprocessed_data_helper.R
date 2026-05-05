rm(list=ls())
graphics.off()

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
files <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9]{8}\\.RData$'))

load(paste0(path_work,tail(files,n=1)))

save(list=c('hemi_left','hemi_right','veri_left','veri_right','desc','Data_comb'),file=paste0(path_work,'workspace_preprocessed_data_',dec,'.RData'))
