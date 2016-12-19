CrossValidation <- function(RES)
{
    require(flexclust)
    require(MASS)
    
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
    rango <- c(0,RES$rango)
    
    FunN <- function(index)
    {Fn <- do.call(fun[index],list(Dist,rango[index]))}
    Fun0 <- function(index)
    {Fn0 <- do.call(fun[index],list(d,rango[index]))}
    
    FuN <- lapply(1:length(fun),FunN)
    names(FuN) <- fun
    
    Mvar <- as.data.frame(matrix(0,ncol=v,nrow=n))
    Mres <- as.data.frame(matrix(0,ncol=v,nrow=n))
    names(Mvar) <- names(RES$data)[3:length(RES$data[1,])]
    names(Mres) <- names(Mvar)
    
    
    for (i in 1:n)
    {
        FuNcv <- FuN
        
        for (j in 1:length(fun))
        {FuNcv[[j]] <- FuNcv[[j]][-i,-i]}
        
        CM <- list()
        
        for (k in 1:length(RES$B))
        {CM[[k]] <- kronecker(FuNcv[[k]],RES$B[[k]])}
        
        names(CM) <- fun
        
        CT <- CM[[1]]
        for (j in 2:length(fun))
        {CT <- CT+CM[[j]]}
        
        #I <- diag(0,v,v)
        #ones <- do.call(rbind, replicate(length(RES$data[,1])-1, I, simplify=FALSE))
        #zeros <- diag(1,v,v)
        #C <- cbind(rbind(CT,t(ones)),rbind(ones,zeros))
        C <- CT
        
        d <- dist2(RES$data[-i,1:2],RES$data[i,1:2])
        
        Fu0 <- lapply(1:length(fun),Fun0)
        
        W <- list()
        c <- list()
        for (j in 1:length(fun))
        {
            c[[j]] <- kronecker(Fu0[[j]],RES$B[[j]])
            #cins <- rbind(c[[j]],I)
            cins <- c[[j]]
            W[[j]] <- ginv(C)%*%cins
        }
        names(W) <- fun
        
        Zor <- (as.vector(t(as.matrix(RES$data[-i,3:length(RES$data[1,])]))))
        Z <- list()
        
        for (j in 1:length(fun))
        {
            zprev2 <- t(W[[j]][1:((n-1)*v),])%*%Zor
            Z[[j]] <- matrix(zprev2,ncol = v,byrow = TRUE)
        }
        names(Z) <- fun
        
        Zt <- Z[[1]]
        varpr <- t(W[[1]][1:((n-1)*v),])%*%c[[1]]
        
        for (j in 2:length(fun))
        {	
            Zt <- Zt+Z[[j]]
            varpr <- varpr+t(W[[j]][1:((n-1)*v),])%*%c[[j]]
        }
        
        colnames(Zt) <- names(RES$data)[3:length(RES$data[1,])]
        
        Mres[i,] <- RES$data[i,3:length(RES$data[1,])]-Zt
        Mvar[i,] <- diag(CT)[1:v]-diag(varpr)
        
    }
    
    CM <- list()
    
    for (k in 1:length(RES$B))
    {CM[[k]] <- kronecker(FuN[[k]],RES$B[[k]])}
    
    names(CM) <- fun
    
    CT <- CM[[1]]
    for (j in 2:length(fun))
    {CT <- CT+CM[[j]]}
    
    C <- CT
    
    d <- dist2(RES$data[,1:2],RES$data[,1:2])
    Fu0 <- lapply(1:length(fun),Fun0)
    
    W <- list()
    c <- list()
    for (j in 1:length(fun))
    {
        c[[j]] <- kronecker(Fu0[[j]],RES$B[[j]])
        #cins <- rbind(c[[j]],I)
        cins <- c[[j]]
        #W[[j]] <- solve(C,cins)
        W[[j]] <- ginv(C)%*%cins
    }
    
    varpr <- t(W[[1]][1:(n*v),])%*%c[[1]]
    
    for (j in 2:length(fun))
    {varpr <- varpr+t(W[[j]][1:(n*v),])%*%c[[j]]}
    
    Var0 <- as.data.frame(t(matrix(diag(C),ncol=n)))
    Var1 <- as.data.frame(t(matrix(diag(varpr),ncol=n)))
    VAR <- colSums(Var0-Var1)
    r2 <- 1-apply(Mres^2,2,mean)/VAR
    rmse <- apply(Mres^2,2,mean)
    
    result <- list(Vartot=VAR,resid=Mres,varres=Mvar,R2=r2,RMSE=rmse)
}