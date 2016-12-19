crossvariogram <- function(coord,values,n=15)
{
require(sp)
data <- cbind(coord,values)
nom <- names(values)
p <- ncol(values)
c <- p*(p+1)/2
comb <- t(cbind(rbind(1:p,1:p),combn(1:p, 2)))

#Calcular variograma experimental
values <- as.matrix(values)
Dist <- as.matrix(dist(coord))
hmin <- Dist[which(Dist == max(Dist), arr.ind = TRUE)][1]/(2.4*n)
h <- seq(1,n)*hmin
delth <- (h[2]-h[1])/4
#h <- c(delth,h)
M0 <- matrix(0,length(values[,1]),length(values[,1]))
M <- list(M0)
eta0 <- rep(0,length(h))
eta <- list(diag(eta0))
N0 <- 0
A0 <- matrix(0,length(values[,1]),length(values[,1]))
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
t <- rep(0,length(h))

for (i in 1:length(h))
{
  t[i] <- t(values[,1])%*%A[[i]]%*%values[,1]
}
Variog <- cbind.data.frame(h,t,unlist(N))
names(Variog) <- c('h','gamma','pairs')
Variog <- list(Variog)

for (j in 2:c)
{
	for (i in 1:length(h))
		{
  		t[i] <- t(values[,comb[j,1]])%*%A[[i]]%*%values[,comb[j,2]]
		}
		Variog[[j]] <- cbind.data.frame(h,t,unlist(N))
		names(Variog[[j]]) <- c('h','gamma','pairs')
}
#Si hay algún Na lo eliminamos
#SUP <- function(t) {Variog[[t]] <- Variog[[t]][-which(is.na(Variog[[1]][,2])),]}
#Variog2 <- lapply(1:c,SUP)
#names(Variog2) <- names(Variog)
#Variog <- Variog2

Variog[[length(Variog)+1]] <- p
Variog[[length(Variog)+1]] <- comb
Variog[[length(Variog)+1]] <- data
names(Variog) <- c(paste(nom[comb[,1]], nom[comb[,2]], sep = " "),"variables","comb","data")
Variog
}