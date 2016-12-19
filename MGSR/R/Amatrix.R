Amatrix <- function(coord,n=15)
{
#Calcular variograma experimental
Dist <- as.matrix(dist(coord))
hmin <- Dist[which(Dist == max(Dist), arr.ind = TRUE)][1]/(2.4*n)
h <- seq(1,n)*hmin
delth <- (h[2]-h[1])/4
M0 <- matrix(0,length(coord[,1]),length(coord[,1]))
M <- list(M0)
eta0 <- rep(0,length(h))
eta <- list(diag(eta0))
N0 <- 0
A0 <- matrix(0,length(coord[,1]),length(coord[,1]))
A <- list(A0)
N <- list(N0)
for (i in 1:length(h)){
  M[[i]] <- M0
  M[[i]][which(Dist >=(h[i]-delth) &  Dist <=(h[i]+delth), arr.ind = TRUE)] <- 1
  N[[i]] <- sum(apply(M[[i]],2,sum))
  eta[[i]] <- diag(eta0)
  eta[[i]] <- diag(apply(M[[i]],2,sum))
  A[[i]] <- (eta[[i]]-M[[i]])/(2*N[[i]])
}
A
}