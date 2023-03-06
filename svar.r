# Stationary Test
Resource<-read.csv(file="E:/숭실대/2019-2/기업경제예측론/GDP_def.csv", col.names=c("GDP", "PRICE"),header=T)
TS<-ts(Resource, start=c(1989,1), end=c(2018,4), frequency=4)
plot(TS)

#Yt(real GDP growth rate)
GDP<- ts(Resource$GDP, start=c(1989,1), end=c(2018,4), frequency=4) 

#Pt(Price change rate)
PRICE<- ts(Resource$PRICE, start=c(1989,1), end=c(2018,4), frequency=4)

#Test H0 : γ(Dickey Fuller)=0 vs H1 : γ≠0
adf.test(GDP, alternative = "stationary")
adf.test(PRICE, alternative = "stationary")



# 2. Reduced form VAR model
VTS<-VAR(TS,type="const", ic="AIC", lag.max=10)
summary(VTS)

# Coefficients matrix of reduced form VAR
A0<-matrix(c(0.62545200,-0.18627009),2,1)
A1<-matrix(c(0.32539609, 0.11211905,-0.16527803, 0.04098555),2,2) > A2<-matrix(c(0.15812592,0.17320110,0.07608662,0.04296556),2,2)
A3<-matrix(c(-0.07909229, 0.24343885, 0.19299639, 0.03019500),2,2)
A4<-matrix(c(-0.02135947, -0.05604426, 0.03547106, 0.34759087),2,2)

TR<-t(Resource)
R<- t(resid(VTS))
E<-matrix(R,2) # Residuals of reduced form VAR
k<-2 #number of endogenous variables
p<-4 #number of lags
t<-120 #number of observations

# Reduce form VAR model
R_VAR<-matrix(0,k,t)
R_VAR[,1]<-TR[,1]
R_VAR[,2]<-TR[,2]
R_VAR[,3]<-TR[,3]
R_VAR[,4]<-TR[,4]
for (i in (p+1) : t) {
+ R_VAR[,i] <-A0 + A1%*%R_VAR[,i-1] + A2%*%R_VAR[,i-2] + A3%*%R_VAR[,i-3] + A4%*%R_VAR[,i-4] + E[,i-p] }



# Cholesky decomposition
# Variance-Covariance matrix
x<-matrix(c(1.4713, -0.2588, -0.2588, 1.0072),2,2)
Ainv<-chol(x)
A<-solve(Ainv)
A

# Structural VAR model
> SV<- matrix(0,k,t)
> SV[,1]<-TR[,1]
> SV[,2]<-TR[,2]
> SV[,3]<-TR[,3]
> SV[,4]<-TR[,4]
> for (i in (p+1) : t) {
+ SV[2,i] <- ((A%*%A0)[2,] +
(A%*%A2%*%SV)[2,i-2] + (A%*%A3%*%SV)[2,i-3] + (A%*%A4%*%SV)[2,i-4] + (A%*%E)[2,i-p])/A[2,2]
+ SV[1,i] <- ((A%*%A0)[1,] - A[1,2]*SV[2,i] + (A%*%A1%*%SV)[1,i-1] + (A%*%A2%*%SV)[1,i-2] + (A%*%A3%*%SV)[1,i-3] + (A%*%A4%*%SV)[1,i-4] +
(A%*%A1%*%SV)[2,i-1] +
(A%*%E)[1,i-p])/A[1,1] }

# Residuals of structural VAR
e <- matrix(0,2,116)
e[1,] <- (A%*%E)[1,]/A[1,1]
e[2,] <- (A%*%E)[2,]/A[2,2]



# Impulse Response Function
amat<-matrix(NA,2,2)
amat[1,2]<-0
SVAR(VTS, Amat=amat)
phi<-irf(SVAR(VTS, Amat = amat),n.ahead=30)
plot(phi)
