
mu=c(.3, .2, .1);
m=length(mu)
covMat=cbind(c(1, .2, -.3), c(.2, 1, .1), c(-.3, .2, 1))
sig=c(.2, .1, .3)

for(i in c(1:m)){
	for(j in c(1:m)){
		covMat[i, j]=covMat[i, j]*sig[i]*sig[j]
	}
}

q=1/10000
z0=cbind(c(1.1, .9, .7))
sigma=t(chol(covMat))

t=1
timeSteps=1000
dt=t/timeSteps
mu=diag(mu)
z=matrix(ncol=timeSteps, nrow=m)
z[, 1]=z0
n=5 #number of loans
w=matrix(nrow=n, ncol=m)
w[1, 1]=.17
w[1, 2]=.44
w[1,3]=.39
w[2, 1]=.04
w[2, 2]=.47
w[2, 3]=.49
w[3, 1]=.02
w[3, 2]=.87
w[3, 3]=.11
w[4, 1]=.41
w[4, 2]=.19
w[4, 3]=.4
w[5, 1]=.52
w[5, 2]=.31
w[5, 3]=.17


#for(i in c(1:n)){
#	wt=runif(m, 0, 1)
#	wt=wt/sum(wt)
#	wt=round(wt, 2)
#	w[i, ]=wt
#}
#l=round(rgamma(n, 10, 1/100))
#r=round(runif(n, 0, 1), 2)
p=c(1:n)*.005
r=c(.13, .15, .18, .14, .78)


l=c(987, 2104, 1264, 576, 377)
#lambda0=0
lambda0=sum(r*l)
r=c(1:n)*0


el=function(z){
	muV=diag(mu)
	myEl=0
	for(j in c(1:m)){
		myEl=myEl+(z[j]/muV[j])*(1-exp(-muV[j]*t))*(z0[j ]-1)+z[j]*t
	}
	return(myEl)
	#return(t+sum((z/muV)*(1-exp(-muV*t))*(t(z0)-1)))
}
vl=function(z1, z2){
	myVar=0
	myA=diag(mu)
	for(i in c(1:m)){
		for(j in c(1:m)){
			myT=t-(1-exp(-myA[i]*t))/myA[i]-(1-exp(-myA[j]*t))/myA[j]
			myT=myT+(1-exp(-(myA[i]+myA[j])*t))/(myA[i]+myA[j])
			myVar=myVar+myT*((covMat[i, j]*z1[i]*z2[j])/(myA[i]*myA[j]))
		}
	
	}
	return(myVar)
}
ek=rbind((p*l))%*%w
ek2=rbind((p*l*l))%*%w
lambda=sum(r*l)

ELnL=el(ek)*(1+q*(lambda+lambda0))

VLnL=(vl(ek, ek)+el(ek2))*(1+q*(lambda+lambda0))^2+el(ek)*(lambda+lambda0)*(lambda+lambda0)*q

marginalE=function(k){
	marg=p[k]*l[k]*el(w[k, ])*(1+q*lambda0)+el(ek)*q*r[k]*l[k]
	return(marg)

}
marginalV=function(k){
	v=(vl(ek, ek)+el(ek2))
	ex=el(ek)
	marg=p[k]*l[k]*l[k]*el(w[k, ])+p[k]*l[k]*vl(w[k,], ek)
	marg=marg*(1+q*lambda0)^2+p[k]*l[k]*q*lambda0*lambda0*el(w[k, ])
	marg=marg+2*r[k]*l[k]*q*v+r[k]*l[k]*q*q*(lambda+2*lambda0)*v
	marg=marg+r[k]*l[k]*(lambda+2*lambda0)*q*ex
	return(marg/sqrt(VLnL))
}
ex=0
ev=0
for(i in c(1:n)){
	print(marginalV(i)+marginalE(i))
	ex=ex+marginalE(i)
	ev=ev+marginalE(i)+marginalV(i)
}
ELnL+sqrt(VLnL)
ev
ex

vl(c(1, 1, 1), c(1, 1, 1))

for(i in c(1:(timeSteps-1))){
	
	x=cbind(rnorm(m, 0, sqrt(dt)))
	z[, i+1]=z[ ,i]+mu%*%(cbind(c(1:m)/c(1:m))-z[, i])*dt+sigma%*%x
} 
Y=sapply(c(1:m), function(i){cumsum(z[i, ]*dt)})
iY=sapply(c(1:m), function(i){z[i, ]*dt})
plot(c(1:timeSteps)*dt, z[3, ], type='l', ylim=c(.6, 1.3), ylab="Value", xlab="Time")
lines(c(1:timeSteps)*dt, z[2, ], lty=2)
lines(c(1:timeSteps)*dt, z[1, ], lty=3)


plot(c(1:timeSteps)*dt, Y[, 3], type='l', ylim=c(0, 1.3), ylab="Value", xlab="Time")
lines(c(1:timeSteps)*dt, Y[, 2 ], lty=2)
lines(c(1:timeSteps)*dt, Y[ ,1], lty=3)


weightedY=w%*%t(Y)
totalP=sapply(c(1:timeSteps), function(i){weightedY[, i]*p})

plot(c(1:timeSteps)*dt, totalP[1, ],  type='l', ylim=c(0, .03), ylab="Value", xlab="Time")
lines(c(1:timeSteps)*dt, totalP[2, ], lty=2)
lines(c(1:timeSteps)*dt, totalP[3, ], lty=3)
lines(c(1:timeSteps)*dt, totalP[4, ], lty=4)
lines(c(1:timeSteps)*dt, totalP[5, ], lty=5)



ts.plot(t(z), col=c("blue", "red", "green"))
ts.plot(Y, col=c("blue", "red", "green"))
ts.plot(iY, col=c("blue", "red", "green"))