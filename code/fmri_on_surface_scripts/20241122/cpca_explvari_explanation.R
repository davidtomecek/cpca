rm(list=ls(all.names=T))
graphics.off()

# creating a random data matrix
nr <- 128
nc <- 32

M <- matrix(rnorm(nr*nc),nrow=nr,ncol=nc)

# variance and standard deviation of variables in M
M_var <- apply(M,2,var)
M_sd <- apply(M,2,sd)

# pca
pca <- prcomp(M)

# explained variance per component
pca_explvar <- pca$sdev^2/sum(pca$sdev^2)

# what is the relationship between the column variance in data matrix and variance or standard deviation in pca components

sum(pca$sdev^2)
sum(M_sd^2)
abs(sum(pca$sdev^2) - sum(M_sd^2))


sum(apply(pca$x,2,var))
abs(sum(pca$sdev^2) - sum(apply(pca$x,2,var)))

head(apply(pca$x,2,var)/sum(apply(pca$x,2,var)))


# so the total variance of the variables in the data matrix is conserved in the pca rotation (duh)
# but does it also hold in the case of cPCA

library(rsvd)
library(complexlm)

path_work <- '/home/anyzjiri/NUDZ/fmri_on_surface/data/workspaces/20241122/'

dec <- 'fsaverage5'

n_pc_ret <- 5

load(paste0(path_work,'ESO_C03419_20230726_1015_1/',dec,'_cpca_data_helper.RData'))

D <- rbind(data_hilbert_left[veri_left,],data_hilbert_right[veri_right,])
D_var <- apply(D,2,var)
# D_sd <- apply(D,2,sd)

cpca <- rpca(D,k=n_pc_ret,center=F,scale=F)

# in the cpca result, the explained variance corresponds to cpca$sdev^2/cpca$var
cpca$sdev^2/cpca$var

# where the total variance is as in the previous case sum of variances of variables in the data matrix,
# the function var from the complexlm package seems to do the same (see variable D_var)
# let's try to find out how it's done (https://en.wikipedia.org/wiki/Complex_random_variable)

# complex expectation E(z) = E(Real(z)) + i*E(imag(z))
# complex variance is Var(z) = E(Mod(z - E(z))^2) = E(Mod(z)^2) - Mod(E(z))^2

# so, for a variable i from the data matrix D, the variance should be 
i <- 456
sum(Mod(D[,i] - mean(D[,i]))^2)/(length(D[,i])-1)
var(D[,i])

compvar <- function(z){
  sum(Mod(z - mean(z))^2)/(length(z)-1)
}

D_compvar <- apply(D,2,compvar)

hist(D_var - D_compvar,100)
mean(D_var - D_compvar)

# it seems that the formula E(Mod(z - E(z))^2) is correct,
# but in computing the outer expectation, the denominator is n-1, because we already estimated one parameter inside
# it's the same as in the ordinary variance


# now the question is how to estimate the explained variance in the group cPCA