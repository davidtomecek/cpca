rm(list=ls(all.names=T))
graphics.off()

Sys.setenv(RGL_USE_NULL=T)

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

visualize_on_surface_pval <- function(x,surf,color_code,file){

  if(length(x) != nrow(surf$vertices)){
    stop('Length of the p-values must match the number of vertices')
  }

  color_code <- color_code[order(color_code$value,decreasing=T),]

  cols <- rep('gray66',times=length(x))
  for(i in 1:nrow(color_code)){
    cols[x < color_code$value[i] & !is.na(x)] <- color_code$color[i]
  }
  cols[is.na(x)] <- 'gray33'

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

color_code <- data.frame(value=c(0.05,0.01,0.001),
                         color=c('orange','yellow','green'))

path_surf <- paste0('/home/anyzjiri/NUDZ/fmri_on_surface/data/20241122/',dec,'/surf/')

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'
file_work <- list.files(path_work,pattern=paste0('^workspace_cpca4set_',dec,'_test_'))

path_rend <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20241122/'
foldname_rend <- 'cpca4set/tests'
if(!dir.exists(paste0(path_rend,'/',foldname_rend,'/'))){dir.create(paste0(path_rend,'/',foldname_rend,'/'),recursive=T)}

library(freesurferformats)
inflated_left  <- read.fs.surface(paste0(path_surf,'lh.inflated'))
inflated_right <- read.fs.surface(paste0(path_surf,'rh.inflated'))

for(comparison in 1:length(file_work)){

  # loading the precomputed group cPCA p-values
  load(paste0(path_work,'/',file_work[comparison]))

  ind_left  <- 1:(length(veri_left))
  ind_right <- (1+length(veri_left)):(length(veri_left) + length(veri_right))

  for(test in 1:length(Test)){

    for(p_value in 2:length(P_values[[test]])){

      if(grepl(pattern='^p_value',x=names(P_values[[test]])[p_value])){
        # left hemisphere
        # visualizing the modulus
        helper <- rep(NA,times=nrow(inflated_left$vertices))
        helper[hemi_left[veri_left]] <- P_values[[test]][[p_value]][ind_left]
        visualize_on_surface_pval(x=helper,
                                  surf=inflated_left,
                                  color_code=color_code,
                                  file=paste0(path_rend,'/',foldname_rend,'/',dec,'_left_cpca4set_set_',set,'_test_',P_values[[test]]$test,'_',names(P_values[[test]])[p_value],'.html'))

        # right hemisphere
        # visualizing the modulus
        helper <- rep(NA,times=nrow(inflated_right$vertices))
        helper[hemi_right[veri_right] - max(hemi_left)] <- P_values[[test]][[p_value]][ind_right]
        visualize_on_surface_pval(x=helper,
                                  surf=inflated_right,
                                  color_code=color_code,
                                  file=paste0(path_rend,'/',foldname_rend,'/',dec,'_right_cpca4set_set_',set,'_test_',P_values[[test]]$test,'_',names(P_values[[test]])[p_value],'.html'))

      }

      if(grepl(pattern='^mean_diff',x=names(P_values[[test]])[p_value])){
        # left hemisphere
        # visualizing the modulus
        helper <- rep(NA,times=nrow(inflated_left$vertices))
        helper[hemi_left[veri_left]] <- P_values[[test]][[p_value]][ind_left]
        visualize_on_surface(x=helper,
                             surf=inflated_left,
                             file=paste0(path_rend,'/',foldname_rend,'/',dec,'_left_cpca4set_set_',set,'_test_',P_values[[test]]$test,'_',names(P_values[[test]])[p_value],'.html'))

        # right hemisphere
        # visualizing the modulus
        helper <- rep(NA,times=nrow(inflated_right$vertices))
        helper[hemi_right[veri_right] - max(hemi_left)] <- P_values[[test]][[p_value]][ind_right]
        visualize_on_surface(x=helper,
                             surf=inflated_right,
                             file=paste0(path_rend,'/',foldname_rend,'/',dec,'_right_cpca4set_set_',set,'_test_',P_values[[test]]$test,'_',names(P_values[[test]])[p_value],'.html'))

      }

    }
  }
}

