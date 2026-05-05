rm(list=ls(all.names=T))
graphics.off()

dec <- 'fsaverage6'

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

files <- list.files(path_work,pattern=paste0('^workspace_cpca4set_opt_',dec,'_[0-9]{8}\\.RData'))

load(paste0(path_work,tail(files,1)))

rm('hemi_left','hemi_right','veri_left','veri_right','Sets','CPCA')
gc()

magn01 <- sapply(Data_comb,function(x) { mean(Mod(as.vector(x))) })
magn02 <- sapply(Data_comb,function(x) { median(Mod(as.vector(x))) })
magn03 <- sapply(Data_comb,function(x) { quantile(Mod(as.vector(x)),prob=0.75) })
magn04 <- sapply(Data_comb,function(x) { quantile(Mod(as.vector(x)),prob=0.9) })

magn11 <- sapply(Data_comb,function(x){median(apply(Mod(x),1,mean))})
magn12 <- sapply(Data_comb,function(x){median(apply(Mod(x),1,median))})
magn13 <- sapply(Data_comb,function(x){median(apply(Mod(x),1,quantile,probs=0.75))})
magn14 <- sapply(Data_comb,function(x){median(apply(Mod(x),1,quantile,probs=0.9))})

save(list=c('desc','magn01','magn02','magn03','magn04',
                   'magn11','magn12','magn13','magn14'),
     file=paste0(path_work,'workspace_magnitudes_',dec,'_',gsub(pattern='-',replacement='',x=Sys.Date()),'.RData'))
