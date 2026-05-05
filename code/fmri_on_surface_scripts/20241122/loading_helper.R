
png('/home/anyzjiri/NUDZ/fmri_on_surface/images/loadings.png',width=1.5e3,height=3e2)
plot(x=seq(from=0,by=0.5,length.out=nrow(CPCA[[1]]$rotation)),
     y=Mod(CPCA[[1]]$rotation[,1]),
     xlab='Time [s]',ylab='Loading cPC 1 - modulus',
     type='l',lwd=2,frame=F,las=2)
abline(v=seq(from=0,by=nrow(CPCA[[1]]$rotation)/length(Sets[[1]]),length.out=length(Sets[[1]])),lty=2,lwd=1.5,col='red')
dev.off()
