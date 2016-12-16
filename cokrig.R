cokrig <- function(RES,xygrid)
{
  require(flexclust)

  Sph <- function(D,as){Fsph <- (ifelse(D>as,Fsph <- 0,Fsph <- (1-1.5*(D/as)+0.5*(D/as)^3)))}#MIRAR SI ES 1-1.5...
  Exp <- function(D,as){Fexp <- (exp(-D/as))}
  Gau <- function(D,as){Fgau <- (exp(-(D/as)^2))}
  Lin <- function(D,as){Flin <- ifelse(D>as,Flin <- 0,Flin <- (1-(D/as)))} #lineal
  Pow <- function(D,as){Fpow <- if(as<=2 & as>0){(1-D^as)}else{stop('rango incorrecto')}}
  Nug <- function(D,as){Fnug <- ifelse(D!=0,Fnug <- 0, Fnug <- 1)}

  v <- dim(RES$B[[1]])[1]
  n <- length(RES$data[,1])
  comb <- t(cbind(rbind(1:v,1:v),combn(1:v, 2)))
  vect <- 1:v
  fun <- colnames(RES$FUN)

  Dist <- as.matrix(dist(RES$data[,1:2]))
  d <- dist2(RES$data[,1:2],xygrid)
  rango <- c(0,RES$rango)

  FunN <- function(index)
  {Fn <- do.call(fun[index],list(Dist,rango[index]))}

  Fun0 <- function(index)
  {Fn0 <- do.call(fun[index],list(d,rango[index]))}

  FuN <- lapply(1:length(fun),FunN)
  names(FuN) <- fun
  Fu0 <- lapply(1:length(fun),Fun0)
  names(Fu0) <- fun
  CM <- list()


  for (i in 1:length(RES$B))
  {
    #C[[i]] <- RES$Bvector[i,index]*FuN[[index]]
    CM[[i]] <- kronecker(FuN[[i]],RES$B[[i]])
  }


  names(CM) <- fun


  CT <- CM[[1]]
  for (i in 2:length(fun))
  {CT <- CT+CM[[i]]}

  #Calculo traza para Cross-Validation
  TR <- n-2*sum(diag(CT))


  I <- diag(0,v,v)
  ones <- do.call(rbind, replicate(length(RES$data[,1]), I, simplify=FALSE))
  zeros <- diag(1,v,v)
  C <- cbind(rbind(CT,t(ones)),rbind(ones,zeros))
  W <- list()
  for (i in 1:length(fun))
  {
    #c <- RES$Bvector[index,i]*Fu0[[i]]
    c <- kronecker(Fu0[[i]],RES$B[[i]])
    onesp <- do.call(rbind, replicate(length(xygrid[,1]), I, simplify=FALSE))
    c <- rbind(c,t(onesp))
    W[[i]] <- ginv(C)%*%c
  }
  names(W) <- fun

  Zor <- (as.vector(t(as.matrix(RES$data[,3:length(RES$data[1,])]))))
  Zdiag <- diag(Zor,v*n,v*n)#Mejor como vector
  Z <- list()
  traza <- list()

  for (i in 1:length(fun))
  {	#Zprev <- Zdiag%*%W[[i]][1:(n*v),]
    #zprev2 <- apply(Zprev,2,sum)
    zprev2 <- t(W[[i]][1:(n*v),])%*%Zor
    Z[[i]] <- matrix(zprev2,ncol = v,byrow = TRUE)
  }
  names(Z) <- fun
  Zt <- Z[[1]]

  for (i in 2:length(fun))
  {
    Zt <- Zt+Z[[i]]
  }
  colnames(Zt) <- names(RES$data)[3:length(RES$data[1,])]
  #sigma <- apply(Zt,2,var)
  FINAL <- list(Z=cbind.data.frame(xygrid,Zt),Traza=TR)

  #Swap <- function(Z,ord=1)
  #{
  #	if(ord==1)
  #	{
  #		i <- 1:v
  #		j <- 1:length(fun)
  #	}
  #	else
  #	{
  #		i <- 1:length(fun)
  #		j <- 1:v
  #	}
  #swap <-lapply(j, function(j) lapply(i, function(i) Z[[i]][[j]]))
  #names(swap) <- names(Z[[1]])
  #for (i in 1:length(fun))
  #{names(swap[[i]]) <- names(Z)}
  #swap
  #}

  #Z0 <- Swap(Zfun(WE0))



}
