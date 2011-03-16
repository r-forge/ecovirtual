###### Biogeography functions - ECOVIRTUAL PACKAGE 
############# Specie Area relationship ###############
spp.area=function(c , z){
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z))
	curve(expr = c*x^z , from=1, to=10^10, xlab="Area", ylim=c(1,500),
	 ylab="Species number", font.lab=2, lwd=2, col=2, main=paste("c = ",c,"; z = ",z), log="xy")
	}
#par(mfrow=c(2,2))
#spp.area(c = 1.5 , z = .25)
#spp.area(c = 2.1 , z = .25)
iRain=function(Nspp , chuva , abund , tempo){
	spp=paste("sp.",1:Nspp)
	ilha=numeric()
	riq=numeric()
	for(i in 1:tempo){
		ilha=union(unique(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund))),ilha)
		riq[i]=length(ilha)
		}
	plot(0:tempo,c(0,riq),type="l",lwd=2,bty='n',xlab="time",ylab="number of species",
	 font.lab=2,col=2,ylim=c(0,Nspp),main=c("Island richness",paste("(Nspp=",Nspp," ; rain=",chuva,")")))
	abline(h=Nspp,lty=3)
	return(riq)
	}
# faca um teste:
#iRain(Nspp=10, chuva=10, abund=c(10,10,10,10,10,10,10,10,10,10), tempo=10)
###### tamanho da ilha x spp
iCol=function(areas, Nspp, chuva.total, abund, tempo){
	n.ilhas=length(areas)
	spp=paste("sp.",1:Nspp)
	ilhas=paste("Island",1:n.ilhas)
	chuva=round(chuva.total*areas/sum(areas))
	riq=numeric()
	par(mfrow=c(ceiling(sqrt(n.ilhas)),ceiling(sqrt(n.ilhas))))
	for(i in 1:n.ilhas){
		riq[i]=iRain(Nspp,chuva[i],abund,tempo)[tempo]
		}
	names(riq)=ilhas
	mod=lm(log10(riq)~log10(areas))
	x11()
	plot(areas,riq,log="xy",font.lab=2,pch=16,col=2,bty="l",
		main=paste("N Islands=",n.ilhas,"; N spp=",Nspp,"; time=",tempo),
		xlab="Island area",ylab="Number of species",ylim=c(1,Nspp))
	abline(mod,lty=2)
	cat("c=",mod[[1]][1],"z=",mod[[1]][2],"\n")
	return(riq)
	}
#arquip(areas=c(10,20,40,80),Nspp=1000,chuva.total=100,abund=rep(10,1000),tempo=10)
iColExt=function(Nspp, chuva, abund, tempo, tx.ext){
	spp=paste("sp.",1:Nspp,sep="")
	ilha=numeric()
	riq=numeric()
	ilha=unique(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund)))
	riq[1]=length(ilha)

	for(i in 2:tempo){
		rr=sample(c("V","M"),length(ilha),replace=TRUE,prob=c((1-tx.ext),tx.ext))
		roleta=data.frame(ilha,rr)
		vivos=roleta[roleta$rr=="V",1]
		ilha=union(sample(spp,chuva,replace=TRUE,prob=abund/sum(abund)),vivos)
		riq[i]=length(unique(ilha))
		}
	plot(0:tempo,c(0,riq),type="l",lwd=2,bty='n',xlab="time",ylab="number of species",
	 font.lab=2,col=2,ylim=c(0,Nspp),
	 main=c("Island Richness",paste("Nspp =",Nspp," ; rain =",chuva,"; tx.ext = ",tx.ext)))
	abline(h=Nspp,lty=3)
	return(riq)
	}
#iColExt(Nspp=100, chuva=5, abund=rep(100,100), tempo=100, tx.ext=.1)

#### biogeografia de Ilhas
MW=function(areas , dist , P , a=1, b=-.01, c=1, d=-.01){
  par(mfrow=c(1,2))
  E=a+b*areas
  I=c+d*dist
  S=numeric()
  for(i in 1:length(areas)){S[i]=I[i]*P/(I[i]+E[i])}
  Tn=numeric()
	for(i in 1:length(areas)){Tn[i]=I[i]*E[i]/(I[i]+E[i])}
  
  curve(I[1]-I[1]*x/P,0,P,bty="n",xaxt="n",yaxt="n",xlab="Species number",
        ylab="Taxas",font.lab=2,lwd=2,ylim=c(0,1))
  curve((E[1]/P)*x,0,P,lwd=2,add=TRUE)
  abline(v=0)
  abline(h=0)
  mtext("P",side=1,at=P,font=2)
  mtext("I1",side=2,at=I[1],font=2,las=1)
  mtext("E1",side=4,at=E[1],font=2,las=1)
  mtext("S1",side=1,at=S[1],font=2,col=2)
  mtext("T1",side=2,at=Tn[1],font=2,las=1,col=2)
  points(S[1],Tn[1],col=2,pch=16,cex=1.3)
  segments(S[1],Tn[1],S[1],0,lty=3,col=2)
  segments(S[1],Tn[1],0,Tn[1],lty=3,col=2)
  
  for(i in 2:length(areas)){
    curve(I[i]-I[i]*x/P,0,P,lwd=2,ylim=c(0,1),add=TRUE,lty=i)
    curve((E[i]/P)*x,0,P,lwd=2,add=TRUE,lty=i)
    mtext(paste("I",i,sep=""),side=2,at=I[i],font=2,las=1)
    mtext(paste("E",i,sep=""),side=4,at=E[i],font=2,las=1)
    mtext(paste("S",i,sep=""),side=1,at=S[i],font=2,col=i+1)
    mtext(paste("T",i,sep=""),side=2,at=Tn[i],font=2,las=1,col=i+1)
    points(S[i],Tn[i],col=i+1,pch=16,cex=1.3)
    segments(S[i],Tn[i],S[i],0,lty=3,col=i+1)
    segments(S[i],Tn[i],0,Tn[i],lty=3,col=i+1)
  }
  
  ex=data.frame(areas=areas,spp=S,dist=dist)
  y=lm(S~areas)[[1]][1]
  z=lm(S~areas)[[1]][2]	
  plot(spp~areas,data=ex,log="xy",font.lab=2,pch=as.character(1:length(areas)),col=2,bty="l",
       main=c("Equilibrium",paste("c = ",round(y,2),"; z = ",round(z,2))),
       xlab="Island area",ylab="Species number",ylim=c(1,P))
  abline(lm(log10(spp)~log10(areas),data=ex),lty=3)
  
  return(ex)
  par(mfrow=c(1,2))
}
#################################################
### Island Biog. Plus Rescue Effect and Internal Colonization  
MW.2.0=function(areas , dist , P , peso.A=.5 , a=1, b=-.01, c=1, d=-.01, e=0, f=.01, g=0, h=.01){
	E=((a+b*areas)*peso.A+(g+h*dist)*(1-peso.A))/(peso.A+(1-peso.A))
	I=((c+d*dist)*peso.A+(e+f*areas)*(1-peso.A))/(peso.A+(1-peso.A))
	S=numeric()
	for(i in 1:length(areas)){S[i]=I[i]*P/(I[i]+E[i])}
	Tn=numeric()
	for(i in 1:length(areas)){Tn[i]=I[i]*E[i]/(I[i]+E[i])}

	curve(I[1]-I[1]*x/P,0,P,bty="n",xaxt="n",yaxt="n",xlab="Species number",
	 ylab="Taxas",font.lab=2,lwd=2,ylim=c(0,1),main="Equilibrium")
	curve((E[1]/P)*x,0,P,lwd=2,add=TRUE)
	abline(v=0)
	abline(h=0)
	mtext("P",side=1,at=P,font=2)
	mtext("I1",side=2,at=I[1],font=2,las=1)
	mtext("E1",side=4,at=E[1],font=2,las=1)
	mtext("S1",side=1,at=S[1],font=2,col=2)
	mtext("T1",side=2,at=Tn[1],font=2,las=1,col=2)
	points(S[1],Tn[1],col=2,pch=16,cex=1.3)
	segments(S[1],Tn[1],S[1],0,lty=3,col=2)
	segments(S[1],Tn[1],0,Tn[1],lty=3,col=2)

	for(i in 2:length(areas)){
		curve(I[i]-I[i]*x/P,0,P,lwd=2,ylim=c(0,1),add=TRUE,lty=i)
		curve((E[i]/P)*x,0,P,lwd=2,add=TRUE,lty=i)
		mtext(paste("I",i,sep=""),side=2,at=I[i],font=2,las=1)
		mtext(paste("E",i,sep=""),side=4,at=E[i],font=2,las=1)
		mtext(paste("S",i,sep=""),side=1,at=S[i],font=2,col=i+1)
		mtext(paste("T",i,sep=""),side=2,at=Tn[i],font=2,las=1,col=i+1)
		points(S[i],Tn[i],col=i+1,pch=16,cex=1.3)
		segments(S[i],Tn[i],S[i],0,lty=3,col=i+1)
		segments(S[i],Tn[i],0,Tn[i],lty=3,col=i+1)
		}

	ex=data.frame(areas=areas,spp=S,dist=dist)
	y=lm(S~areas)[[1]][1]
	z=lm(S~areas)[[1]][2]

	plot(spp~areas,data=ex,log="xy",font.lab=2,pch=16,col=2,bty="l",
	 main=c("Species-area relationship",paste("c = ",round(y,2),"; z = ",round(z,2))),
	 xlab="Island area",ylab="Species number",ylim=c(1,P))
	abline(lm(log10(spp)~log10(areas),data=ex),lty=3)

	return(ex)
	}
########### MODELOS NULO ###############
rand.walk <- function(n=1,step=1,ciclo=1e5,cont=1e3,x1=NULL){
  if(is.null(x1)){
    x1 <- sample(1:200,n,replace=TRUE)
  }
  results <- matrix(NA,nrow=1+ciclo/cont,ncol=n) 
  results[1,] <- x1
  X <- x1
  for(i in 2:(1+ciclo/cont)){
    for(j in 1:cont){
      X[X<=0] <- NA
      X <- X +sample(c(step,-1*step),n,replace=TRUE)
    }
    results[i,] <- X
  }
  results[is.na(results)] <- 0
  time <- seq(0,ciclo,by=cont)
  matplot(time,results,type="l", col=rainbow(n),lwd=2, xlab="Steps", ylab="Distancia from the edge")
  abline(h=0,lwd=4, main="Randon Walk")
}
#rand.walk(n=10,step=10,ciclo=1e4,cont=1e3)
#rand.walk(n=10,step=2,ciclo=1e4,cont=1e3)
#rand.walk(n=10,step=2,ciclo=5e4,cont=1e3)
##### Game
ext.game <- function(aposta=1,total=100){
  X <- total/2
  results <- X
  while(X>0&X<total){
    X <- X+sample(c(aposta,-1*aposta),1)
    results <- c(results,X)
  }
  plot(1:length(results),results, type="l", col="blue",ylim=c(0,total), xlab="Cicle", ylab="Number of Individuos")
  lines(1:length(results),total-results, col="red")
  abline(h=c(0,total),lty=2)
  legend(total*0.5,dim(results)*0.8, legend=c("sp 1", "sp2"), lty=1, col=c("red", "blue")  )
}
#old<-par(mfrow=c(2,2))
#ext.game(aposta=1,total=20)
#ext.game(aposta=1,total=50)
#ext.game(aposta=1,total=100)
#ext.game(aposta=1,total=200)
#par(old)
##Modelo neutro sem imigracao 
sim.hub1=function(S= 100, j=10, D=1, ciclo=1e4, step=1000){ 
  ## Tamanho da comunidade
  rich <- function(x)length(unique(x))
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##CONDICOES INICIAIS##
  ##Deduzidas de acordo com o modelo de Hubbell:
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ind.mat[,1] <- rep(1:S,each=j)
  cod.sp <- ind.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    for(j in 1:step){
      ##Indice dos individuos que morrem
      morte <- sample(1:J,D)
      ##Indice dos individuos que produzem filhotes para substituir os mortos
      novos <- sample(1:J,D,replace=TRUE)
      ##Substituindo
      cod.sp[morte]<-cod.sp[novos]
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  invisible(ind.mat)
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (ciclos)", ylab="Number of species",ylim=c(0,S), type="l", main=paste("Neutral Model Without Colonization", "\n S=",S," J=",J), sub=paste("Mean extintion=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"), cex.sub=0.7) 
}

#par(mfrow=c(2,2))
#sim.hub1(j=2,ciclo=2e4,step=1e2)
#sim.hub1(j=5,ciclo=2e4,step=1e2)
#sim.hub1(j=10,ciclo=2e4,step=1e2)
#sim.hub1(j=20,ciclo=2e4,step=1e2)
#par(mfrow=c(1,1))
##
## Com migracao de uma metacomunidade com a composicao inicial
sim.hub2=function(S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01){ 
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step) 
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ## Rotulo de especies para cada um dos inividuos
  ind.mat[,1] <- rep(1:S,each=j)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  cod.sp <- ind.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    for(j in 1:step){
      ##Indice dos individuos que morrem
      morte <- sample(1:J,D)
      ## Indice dos individuos mortos que serao repostos por migrantes
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ##Indice dos individuos que produzem filhotes para substituir os mortos
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:J,sum(defora),replace=TRUE)
      ##Substituindo
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-ind.mat[,1][novosf]
      }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  invisible(ind.mat)
  plot(tempo,apply(ind.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutral Dynamics - Original Community Colonization",sub=paste( "S=",S," J=",J," m=",m,"Mean Extintion rate =",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"),ylim=c(0,S), cex.sub=0.7)
  }
#teste2 <- sim.hub2(j=2,ciclo=2e4,step=1e2,m=0.1)
## Com migracao de uma metacomunidade com especiacao
sim.hub3=function(Sm=200, jm=20, S= 100, j=10, D=1, ciclo=1e4, step=1000, m=0.01, nu=0.001){ 
  ## Tamanho da metacomunidade
  Jm <- Sm*jm
  ## Tamanho da comunidade
  J <- S*j
  ##Matrizes para guardar os resultados
  ## matriz da especie de cada individuo por ciclo
  ## Na metacomunidade
  meta.mat=matrix(nrow=Jm,ncol=1+ciclo/step) 
  ## Na comunidade
  ind.mat=matrix(nrow=J,ncol=1+ciclo/step)
  ##CONDICOES INICIAIS##
  ## Todas as especies comecam com o mesmo numero de individuos (j=J/S)
  ## METACOMUNIDADE
  meta.mat[,1] <- rep(1:Sm,each=jm)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  meta.sp <- meta.mat[,1]
  ##COMUNIDADE
  ## Rotulo de especies para cada um dos individuos
  ind.mat[,1] <- rep(1:S,each=j)
  ## Repetindo este rotulo no vetor que sofrera modificacoes
  cod.sp <- ind.mat[,1]
  ##Aqui comecam as simulacoes
  for(i in 2:(1+ciclo/step)){
    for(j in 1:step){
      ##Indice dos individuos que morrem
      ## Na comunidade
      morte <- sample(1:J,D)
      ## Na metacomunidade
      meta.morte <- sample(1:Jm,D)
      ## Indice dos individuos mortos da comunidade que serao repostos por migrantes 
      defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(m,1-m))
      ## Indice dos individuos mortos da metacomunidade que serao repostos por novas especies 
      meta.defora <- sample(c(TRUE,FALSE),size=D,replace=TRUE,prob=c(nu,1-nu))
      ##Indice dos individuos que produzem filhotes para substituir os mortos da comunidade
      novosd <- sample(1:J,D-sum(defora),replace=TRUE)
      novosf <- sample(1:Jm,sum(defora),replace=TRUE)
      ##Indice dos individuos que produzem filhotes para substituir os mortos da metacomunidade
      meta.novosd <- sample(1:Jm,D-sum(meta.defora),replace=TRUE)
      meta.novosf <- sample(1:Jm,sum(meta.defora),replace=TRUE)
      ##Substituindo
      ## N metacomunidade ##
      ## Mortos por propagulos de dentro
      if(length(meta.novosd)>0){
        meta.sp[meta.morte[!meta.defora]]<-meta.sp[meta.novosd]
      }
      ## Mortos por novas especies
      if(length(meta.novosf)>0){
        meta.sp[meta.morte[meta.defora]]<-max(meta.sp)+1
      }
      ## Na comunidade ##
      ## Mortos por propagulos de dentro
      if(length(novosd)>0){
        cod.sp[morte[!defora]]<-cod.sp[novosd]
      }
      ## Mortos por propagulos de fora
      if(length(novosf)>0){
        cod.sp[morte[defora]]<-meta.sp[novosf]
      }
    }
    ## A cada step ciclos os resultados sao gravados
    ind.mat[,i] <- cod.sp
    meta.mat[,i] <- meta.sp
  }
  tempo <- seq(0,ciclo,by=step)
  colnames(ind.mat) <- tempo
  colnames(meta.mat) <- tempo
  resultados <- list(metacomunidade=meta.mat,comunidade=ind.mat)
  ## Graficos
  plot(tempo,apply(meta.mat,2,rich), xlab="Time (cicles)", ylab="Number of species", type="l",
       main="Neutra Dynamics - Metacomunity Colonization" ,sub=paste( "Jm=",Jm," nu=",nu," Theta=",2*Jm*nu, "S=",S," J=",J," m=",m, " Mean Extintion Rate=",(S-rich(ind.mat[,ncol(ind.mat)]))/ciclo,"sp/ciclo"), col="blue",  ylim=c(0,max(apply(meta.mat,2,rich))), cex.sub=0.7)
  lines(tempo,apply(ind.mat,2,rich),col="red")
  text(ciclo/2 ,length(unique(ind.mat[,round(dim(ind.mat)[2]/2)]))*1.3, "Metacommunity", col="red")
  text(ciclo/2 ,length(unique(meta.mat[,round(dim(ind.mat)[2]/2)]))*0.95, "Community", col="blue")
  invisible(resultados)
  }
#teste3 <- sim.hub3(j=2, ciclo=2e4,step=1e2,m=0.1)
#teste3 <- sim.hub3(j=2, ciclo=2e5,step=1e3,nu=0.00001,m=0.1)

