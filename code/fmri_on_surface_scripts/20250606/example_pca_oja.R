rm(list=ls(all.names=T))
graphics.off()


# number of observations
n <- 128

# number of variables
p <- 32

# random data matrix
D <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
D_sc <- apply(D,2,scale,scale=T,center=T)

# reference pca solution
pca_ref <- prcomp(D,center=T,scale=T)

# an attempt to perform the Oja method

# initializing the first loading vector and normalizing it to unit length
pc1 <- matrix(rnorm(n=p),nrow=p,ncol=1)
pc1 <- pc1/sqrt(sum(pc1^2))

pc1_sha <- pc1
pc1_oja <- pc1

# performing several iterations of the Oja method
niter <- 1024

# learning rate
nu <- 0.0005

for(ni in 1:niter){
  for(no in 1:n){
    # shamir 2016
    pc1_sha <- (diag(p) + nu*(D_sc[no,]%*%t(D_sc[no,]))) %*% pc1_sha
    # pc1_sha <- pc1_sha/sqrt(sum(pc1_sha^2))
    
    # oja & karhunen 1985
    A <- D_sc[no,]%*%t(D_sc[no,])
    pc1_oja <- pc1_oja + nu*(A %*% pc1_oja - as.numeric(t(pc1_oja) %*% A %*% pc1_oja) * pc1_oja)
    # pc1 <- pc1/sqrt(sum(pc1^2))
  }
}

# renormalizing the prinicipal component
pc1_sha <- pc1_sha/sqrt(sum(pc1_sha^2))


# comparing the results
cbind(pca_ref$rotation[,1],pc1_oja,pc1_sha)
cor(cbind(pca_ref$rotation[,1],pc1_oja,pc1_sha))
# and also visually

yupli <- max(abs(c(pca_ref$rotation[,1],pc1_oja,pc1_sha)))

plot(x=1:p,
     y=pca_ref$rotation[,1],
     ylim=c(-yupli,yupli),
     type='l',col='green',lwd=2,
     frame=F)

lines(x=1:p,y= pc1_oja,col='red',lwd=1.2)
lines(x=1:p,y=-pc1_oja,col='red',lwd=1.2,lty=2)

lines(x=1:p,y= pc1_sha,col='blue',lwd=1.2)
lines(x=1:p,y=-pc1_sha,col='blue',lwd=1.2,lty=2)

legend('bottom',legend=c('SVD','Oja & Karhunen 1985','Shamir 2016'),col=c('green','red','blue'),lwd=c(2,1.2,1.2),horiz=T,bty='n')

# estimating more eigenvectors

library(pracma)

# initializing the eigenvectors
r <- 5
R <- matrix(runif(p*r,min=0,max=1),nrow=p,ncol=r)
R <- apply(R,2,scale,scale=T,center=T)

# learning rate and number of iterations
nu <- 0.001
niter <- 1024

# block size
bs <- r^2

for(ni in 1:niter){
  # this seems to work
  # computations based on one observation
  for(no in 1:n){
    A <- D_sc[no,] %*% t(D_sc[no,])
    # R <- R + nu*(A %*% R - (t(R) %*% A %*% R) %*% t(R))

    R <- (diag(p) + nu*A) %*% R

    # # https://meyn.ece.ufl.edu/wp-content/uploads/sites/77/archive/spm_files/Oja08/VTOJA08.pdf
    # # not sure about this
    # a_n <- (1/(1+ni)) * (1/(1+sum(diag(R%*%t(R)))))
    # xi_n <- 0.5*(ni+1)*(ni)*gammaz(0.5*(ni+1))*zeta(ni+1)/sqrt(pi^(ni+1))
    #
    # R <- R + a_n * (diag(p) - R%*%t(R)) %*% A %*% R + a_n*xi_n

    R <- qr.Q(qr(R))
  }

  # # computations based on more observations
  # si <- sample(x=1:nrow(D_sc),size=bs)
  # A <- (t(D_sc[si,]) %*% (D_sc[si,]))/bs
  # R <- (diag(p) + nu*A) %*% R
  # R <- qr.Q(qr(R))
}

# comparing the results

for(i in 1:r){
  print(cbind(pca_ref$rotation[,i],R[,i]))
  print(cor(x=pca_ref$rotation[,i],y=R[,i]))
}

cor(pca_ref$rotation[,1:r],R)

diag(cor(pca_ref$rotation[,1:r],R))



























# next try to work out how to do the same with complex numbers

# first the first principal component case
library(complexlm)

# number of observations
n <- 512

# number of variables
p <- 64

# random data matrix
D <- matrix(complex(real=rnorm(n=n*p),imaginary=rnorm(n=n*p)),nrow=n,ncol=p)
D <- scale(D,center=T,scale=T)

# reference pca solution
pca_ref <- prcomp(D,center=F,scale=F)

# number of components
r <- 3

# # initializing the first loading vector and normalizing it to unit length
# R <- matrix(complex(real=runif(n=p*r,min=0,max=1),imaginary=runif(n=p*r,min=0,max=1)),nrow=p,ncol=r)
# R <- qr.Q(qr(R))

# trying to get a better initialization
# number of samples for initialization
s <- 128
pca_init <- svd(D[sample(1:n,size=s),])
R <- pca_init$v[,1:r]

# keeping the initial values for plotting
R_init <- R

# performing several iterations of the Oja method
niter <- 1e4

# learning rate
nu <- 0.001

# block size
b <- 64

# doing the iterations  
for(ni in 1:niter){
  # for(no in 1:n){
  #   A <- D[no,] %*% Conj(t(D[no,]))
  #   R <- (diag(p) + nu*A) %*% R
  #   R <- qr.Q(qr(R))
  # }
  
  bi <- sample(1:n,size=b)
  A <- (Conj(t(D[bi,])) %*% ((D[bi,])))/b
  R <- (diag(p) + nu*A) %*% R
  R <- qr.Q(qr(R))
  
}

# comparing the results

for(u in 1:r){
  # print(cbind(pca_ref$rotation[,u],R[,u]))
  print(Mod(cor(x=pca_ref$rotation[,u],y=R[,u])))

  yupli <- max(Mod(c(pca_ref$rotation[,u],R[,u])))

  plot(x=1:p,
       y=Mod(pca_ref$rotation[,u]),
       ylim=c(0,yupli),
       type='l',col='green',lwd=2,
       main=paste('PC',u),
       frame=F)

  lines(x=1:p,y=Mod(R[,u]),col='red',lwd=1.2)
  lines(x=1:p,y=Mod(R_init[,u]),col='green',lwd=1.2,lty=2)
  
  legend('bottom',legend=c('SVD','Oja method'),col=c('green','red'),lwd=c(2,1.2),horiz=T,bty='n')
}



# trying the incremental circular principal component analysis by Papaioannou


# number of observations
n <- 512

# number of variables
p <- 64

# random data matrix
D <- matrix(complex(real=rnorm(n=n*p),imaginary=rnorm(n=n*p)),nrow=n,ncol=p)
D <- scale(D,center=T,scale=T)

# reference pca solution by prcomp
pca_ref <- prcomp(D,center=F,scale.=F,retx=F)

# the incremental cPCA requires orthonormal basis and corresponding eigenvalues from previous iteration
# i am going to initialize it by pca on smaller subset of the data

# number of samples used for the initialization
m <- 64

# initialization pca/svd
pca_init <- svd(D[sample(x=1:n,size=m),])

# number of components
r <- 5

# initial values of the projection matrix and the eigenvalues
U <- pca_init$u[1:r,]
l <- pca_init$d[1:r]

# number of samples in the incremental data
o <- 16

for(i in 1){
  
  Xa <- D[sample(1:n,size=o),]
  
  Za <- Xa %*% Conj(t(U))
  
  Ua <- qr.Q(qr(    t(  Za -  Za %*% (((U)) %*% Conj(t(U)))  )    ),ncol=m)
  
  R <- cbind(rbind(diag(l),
                   matrix(0,nrow=2*o-r,ncol=r)),
             rbind(Conj(t(Conj(t(U)) %*% Conj(t(Xa)))),
                   Conj(t(Conj(t(Ua)) %*% Conj(t(Xa))))))
  
  Rsvd <- svd(R)
  
  U <- cbind(U,Ua) %*% Rsvd$u
  l <- Rsvd$d[1:r]
  
}

