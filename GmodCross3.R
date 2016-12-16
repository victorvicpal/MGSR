GmodCross3 <- function(CV,fun,a)
{
	h <- CV[[1]][,1]
	m <- length(fun)
	cfun <- rep(0,m)
	b <- rep(0,m)
	nug <- rep(0,m)
	pin <- CV$variables*(1+CV$variables)/2


	Sph <- function(a,h){
    Fsph <- ifelse(h>a,Fsph <- 1,Fsph <- (1.5*(h/a)-0.5*(h/a)^3))
    Fsph
	}

	Exp <- function(a,h){
    Fexp <- (1-exp(-h/a))
    Fexp
	}

	Gau <- function(a,h){
  	Fgau <- (1-exp(-(h/a)^2))
  	Fgau
	}

	Lin <- function(a,h){
	Flin <- ifelse(h>a,Flin <- 1,Flin <- h/a)
	Flin
	}

	Pow <- function(a,h){
	Fpow <- h^a
	Fpow
	}

	FUN <- list(0)
	for (i in 1:length(fun))
	{
		if (fun[i]=='Exp')
		{FUN[[i]]<-Exp(a[i],h)}
		else {if (fun[i]=='Sph'){FUN[[i]]<-Sph(a[i],h)}
			  else {if (fun[i]=='Gau'){FUN[[i]]<-Gau(a[i],h)} 
			  else {if (fun[i]=='Lin'){FUN[[i]]<-Lin(a[i],h)}
			  else {if (fun[i]=='Pow'){FUN[[i]]<-Pow(a[i],h)} else{stop("Función incorrecta")}}}}}}
			  
	names(FUN) <- fun


	ae <- rep(0,length(which(fun=='Exp')))
	Exp <- ae
	ie <- 1
	ag <- rep(0,length(which(fun=='Gau')))
	Gau <- ag
	ig <- 1
	as <- rep(0,length(which(fun=='Sph')))
	Sph <- as
	is <- 1
	al <- rep(0,length(which(fun=='Lin')))
	Lin <- al
	il <- 1
	ap <- rep(0,length(which(fun=='Pow')))
	Pow <- ap
	ip <- 1
	#Rangos
	#Ojo! Diferentes rangos para cada función
	for(index in 1:m)
	{
		if (fun[index]=='Exp')
		{ae[ie] <- a[index]
			ie <- ie+1}
		else{if(fun[index]=='Gau'){ag <- a[index]
									ig <- ig+1}
			else{if(fun[index]=='Sph'){as[is] <- a[index]
										is <- is+1}
				else{if(fun[index]=='Pow'){ap[ip] <- a[index]
											ip <- ip+1}
					else{al[il] <- a[index]
						il <- il+1}
					}
				}
			}
	}

	#Función/Funciones
	for(i in 1:length(ae)){Exp[i] <- paste("I(exp(-h/",ae[i],"))",sep="")} #Cortar los CV para a>h
	for(i in 1:length(ag)){Gau[i] <- paste("I(exp(-(h/",ag[i],")^2))",sep="")}
	for(i in 1:length(as)) {Sph[i] <- paste("I((1.5*(h/",as[i],")-0.5*(h/",as[i],")^3))",sep="")}
	for(i in 1:length(al)){Lin[i] <- paste("I(h/",al[i],")",sep="")} #Cortar los CV para a>h
	for(i in 1:length(ap)){Pow[i] <- paste("I(h^",ap[i],")",sep="")}

	ie <- ifelse(ie==1,ie <- 0, ie <-1)
	ig <- ifelse(ig==1,ig <- 0, ig <-1)
	is <- ifelse(is==1,is <- 0, is <-1)
	ip <- ifelse(ip==1,ip <- 0, ip <-1)
	il <- ifelse(il==1,il <- 0, il <-1)

	for(index in 1:m)
	{
		if (fun[index]=='Exp')	{cfun[index] <- Exp[ie]
								ie <- ie+1}
		else{if(fun[index]=='Gau')	{cfun[index] <- Gau[ig]
								ig <- ig+1}
			else{if(fun[index]=='Sph'){cfun[index] <- Sph[is]
										is <- is+1}
				else{if(fun[index]=='Pow'){cfun[index] <- Pow[ip]
											ip <- ip+1}
					else{if(fun[index]=='Lin'){cfun[index] <- Lin[il]
												il <- il+1}}
					}
				}
			}
	}


	if(m==1){formula <- as.formula(paste("gamma~",cfun))
	}else{formula <- as.formula(paste("gamma~",paste(cfun,collapse = '+')))}

	cvfit <- function(index)
	{fit <- if((Reduce("|",fun=='Sph')|Reduce("|",fun=='Lin'))==TRUE)
		{
			L <- max(al,as)
			fit <- glm(formula,data=CV[[index]][h<L,1:2])
		}
		else{fit <- glm(formula,data=CV[[index]][,1:2])}
	fit$coefficients
	}

	FIT <- lapply(1:pin,cvfit)
	names(FIT) <- names(CV)[1:pin]

	Bvector <- matrix(0,length(FIT),length(FIT[[1]]))
	Mfun <- matrix(0,)

	for (i in 1:length(Bvector[,1]))
	{Bvector[i,] <- FIT[[i]]}

	colnames(Bvector) <- c('nug',fun)
	rownames(Bvector) <- names(FIT)

	r <- rep(1,length(FUN[[1]]))

	for (i in 1:length(fun))
	{r <- cbind(r,unlist(FUN[[i]]))}
	colnames(r) <- c('Nug',fun)

	RESULT <- list(Bvector=Bvector,FUN=r)

}