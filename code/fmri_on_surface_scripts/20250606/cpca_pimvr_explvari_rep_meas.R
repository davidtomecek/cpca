rm(list=ls(all.names=T))
graphics.off()

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20250606/'
load(paste0(path_work,'workspace_fsaverage6_cpca_pimvr_partial_explvari.RData'))

path_imag <- '/home/anyzjiri/NUDZ/fmri_on_surface/images/20250606/explvari/'
if(!dir.exists(path_imag)){dir.create(path_imag,recursive=T)}

for(pc in 1:3){

  png(filename=paste0(path_imag,'repmeas_explvari_cpc_',pc,'.png'),
      width=5,height=7.5,res=300,units='in')
  layout(mat=matrix(c(1,1,2,3),ncol=2,byrow=T),widths=c(1/2,1/2),heights=c(2/3,1/3))
  
  cols <- factor(c(ExplVari[[1]]$group[ExplVari[[1]]$visit=='1'],ExplVari[[1]]$group[ExplVari[[1]]$visit=='2']),levels=c('C','P'))
  levels(cols) <- c('darkblue','firebrick')
  cols <- as.character(cols)
  
  cc_all <- round(cor(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1'],
                      y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2'],
                      method='pearson'),
                  digits=3)
  ulim_all <- ceiling(10*max(ExplVari[[1]][,pc]))/10
  plot(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1'],
       y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2'],
       main=paste0('cPC ',pc,' - All data'),sub=bquote(rho == .(cc_all)),xlab='EV V1',ylab='EV V2',
       frame=F,xlim=c(0,ulim_all),ylim=c(0,ulim_all),pch=19,col=cols)
  legend('bottomright',legend=c('Healthy controls','Patients'),col=c('darkblue','firebrick'),pch=19,box.lty=0,bg='gray95')
  abline(a= 0.00,b=1,lty=2)
  abline(a=-0.05,b=1,lty=2,lwd=0.5)
  abline(a= 0.05,b=1,lty=2,lwd=0.5)
  abline(a=-0.10,b=1,lty=2,lwd=0.3)
  abline(a= 0.10,b=1,lty=2,lwd=0.3)
  
  cc_hc <- round(cor(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='C'],
                     y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='C'],
                     method='pearson'),
                 digits=3)
  ulim_hc <- ceiling(10*max((ExplVari[[1]][,pc])[ExplVari[[1]]$group=='C']))/10
  plot(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='C'],
       y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='C'],
       main=paste0('cPC ',pc,' - Healthy controls'),sub=bquote(rho == .(cc_hc)),xlab='EV V1',ylab='EV V2',
       frame=F,xlim=c(0,ulim_hc),ylim=c(0,ulim_hc),pch=19,cex=0.6,col='darkblue')
  abline(a= 0.00,b=1,lty=2)
  abline(a=-0.05,b=1,lty=2,lwd=0.5)
  abline(a= 0.05,b=1,lty=2,lwd=0.5)
  abline(a=-0.10,b=1,lty=2,lwd=0.3)
  abline(a= 0.10,b=1,lty=2,lwd=0.3)
  
  cc_p <- round(cor(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='P'],
                    y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='P'],
                    method='pearson'),
                digits=3)
  ulim_p <- ceiling(10*max((ExplVari[[1]][,pc])[ExplVari[[1]]$group=='P']))/10
  plot(x=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='P'],
       y=(ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='P'],
       main=paste0('cPC ',pc,' - Patients'),sub=bquote(rho == .(cc_p)),xlab='EV V1',ylab='EV V2',
       frame=F,xlim=c(0,ulim_p),ylim=c(0,ulim_p),pch=19,cex=0.6,col='firebrick')
  abline(a= 0.00,b=1,lty=2)
  abline(a=-0.05,b=1,lty=2,lwd=0.5)
  abline(a= 0.05,b=1,lty=2,lwd=0.5)
  abline(a=-0.10,b=1,lty=2,lwd=0.3)
  abline(a= 0.10,b=1,lty=2,lwd=0.3)
  dev.off()
  
  diff_all <- abs((ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1'] - (ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2'])
  diff_hc  <- abs((ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='C'] - (ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='C'])
  diff_p   <- abs((ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='1' & ExplVari[[1]]$group=='P'] - (ExplVari[[1]][,pc])[ExplVari[[1]]$visit=='2' & ExplVari[[1]]$group=='P'])
  
  png(filename=paste0(path_imag,'repmeas_explvari_hist_cpc_',pc,'.png'),
      width=5,height=7.5,res=300,units='in')
  layout(mat=matrix(c(1,1,2,3),ncol=2,byrow=T),widths=c(1/2,1/2),heights=c(2/3,1/3))
  
  ulim_hist <- ceiling(10*max(diff_all))/10
  hist_all <- hist(x=diff_all,
                   breaks=seq(from=0,to=ulim_hist,by=0.01),
                   main=paste0('cPC ',pc,' - All data'),xlab=expression(EV ~ abs(V2 - V1)),
                   col='gray60')
  hist_hc  <- hist(x=diff_hc,
                   breaks=seq(from=0,to=ulim_hist,by=0.01),
                   main=paste0('cPC ',pc,' - Healthy controls'),xlab=expression(EV ~ abs(V2 - V1)),
                   ylim=c(0,max(hist_all$counts)),col=adjustcolor('darkblue',alpha.f=0.6))
  hist_p   <- hist(x=diff_p,
                   breaks=seq(from=0,to=ulim_hist,by=0.01),
                   main=paste0('cPC ',pc,' - Patients'),xlab=expression(EV ~ abs(V2 - V1)),
                   ylim=c(0,max(hist_all$counts)),col=adjustcolor('firebrick',alpha.f=0.6))
  
  dev.off()
}