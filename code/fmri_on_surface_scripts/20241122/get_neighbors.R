rm(list=ls())
graphics.off()

library(freesurferformats)

dec <- 'fsaverage6'

path_data <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')
path_save <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/',dec,'/')
if(!dir.exists(path_save)){dir.create(path_save,recursive=T)}

# left hemisphere
# surface data
surf_left  <- read.fs.surface(paste0(path_data,'lh.pial'))

# extracting the list of neighbors for each vertex
neigh_left <- vector(length=nrow(surf_left$vertices),mode='list')

for(i in 1:nrow(surf_left$faces)){
  neigh_left[[surf_left$faces[i,1]]] <- c(neigh_left[[surf_left$faces[i,1]]],surf_left$faces[i,c(2,3)])
  neigh_left[[surf_left$faces[i,2]]] <- c(neigh_left[[surf_left$faces[i,2]]],surf_left$faces[i,c(1,3)])
  neigh_left[[surf_left$faces[i,3]]] <- c(neigh_left[[surf_left$faces[i,3]]],surf_left$faces[i,c(1,2)])
}

neigh_left <- lapply(neigh_left,function(x) sort(unique(x)))

save(neigh_left,file=paste0(path_save,'neigh_left.RData'))

rm('surf_left','neigh_left')


# right hemisphere
# surface data
surf_right <- read.fs.surface(paste0(path_data,'rh.pial'))

# extracting the list of neighbors for each vertex
neigh_right <- vector(length=nrow(surf_right$vertices),mode='list')

for(i in 1:nrow(surf_right$faces)){
  neigh_right[[surf_right$faces[i,1]]] <- c(neigh_right[[surf_right$faces[i,1]]],surf_right$faces[i,c(2,3)])
  neigh_right[[surf_right$faces[i,2]]] <- c(neigh_right[[surf_right$faces[i,2]]],surf_right$faces[i,c(1,3)])
  neigh_right[[surf_right$faces[i,3]]] <- c(neigh_right[[surf_right$faces[i,3]]],surf_right$faces[i,c(1,2)])
}

neigh_right <- lapply(neigh_right,function(x) sort(unique(x)))

save(neigh_right,file=paste0(path_save,'neigh_right.RData'))

rm('surf_right','neigh_right')
