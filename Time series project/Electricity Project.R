# Sourse: https://www.kaggle.com/shenba/time-series-datasets?select=Electric_Production.csv
data=read.csv("Electric_Production.csv")
colnames(data)=c("date","unit")
attach(data)
plot(unit,type="l",xlab="Time",main="Electricity Production")

##  A) Testing for the presence of trend using Relative ordering test:
n=length(unit)
count=0
for(i in 1:(n-1)){
for(j in (i+1):n){
if(unit[i]>unit[j]){
count=count+1
}}}
count
tau1=1-((4*count)/(n*(n-1)))
tau1
var_tau1=(2*(2*n+5)/(9*n*(n-1)))
z=tau1/sqrt(var_tau1)
z
qnorm(1-0.025)
#  Reject H0 and conclude that trend is present

## B) Test for present of Seasonality Friedman Test

Yt=unit
t=1:length(Yt)
Xt=Yt-predict(lm(Yt~t))
plot(Xt,type="l",xlab="Time",main="Electricity Consumption")
Xt=Xt[-397]
length(Xt)
A=matrix(Xt,nrow=12,byrow=F)
Mi=apply(A,1,sum)
c=33
r=12
X=12/(c*r*(r+1))*sum((Mi-c*(r+1)/2)^2)
X
qchisq(1-0.05,r-1)
# H0 is rejected and conclude that seasonality is present

## C) Determination of Period of Seasonality
D12.Xt=Xt[13:396]-Xt[1:384]
plot(D12.Xt,type="l",xlab="Time",main="Electricity Consumption")

u=0
rpm=D12.Xt
for(i in 2:(384-1))
{
if(((rpm[i]>rpm[i+1])&&(rpm[i]>rpm[i-1]))||((rpm[i]<rpm[i+1])&&(rpm[i]<rpm[i-1])))
{
u=u+1
}
}
u
e_u=2*(n-2)/3
v_u=(16*n-29)/90
z=(u-e_u)/v_u
z
qnorm(1-0.025)

# Since D12 makes he series purely random so the period of seasonalaty is 12

### So, Let us move to estimation of components of time series

#Rough Estimation of trend using 12 point moving average
Yt=Yt[-397]
Y1t=array(0)
for(i in (1:385)){
	Y1t[i]=mean(Yt[i:(i+11)])
}

j=1
a1=array(0)
for(i in 7:2){
	a1[j]=(i*Yt[1]+sum(Yt[2:(12-(i-1))]))/12
	j=j+1
}
length(a1)

j=1
a2=array(0)
for(i in 2:7){
	a2[j]=(i*Yt[396]+sum(Yt[395:(385+(i-1))]))/12
	j=j+1
}
length(a2)

Y11t=c(a1,Y1t,a2)
length(Y11t)

m1t=array(0)
for(i in (1:396)){
	m1t[i]=mean(Y11t[i:(i+1)])
}

length(m1t)   # m1t is the moving average values, which are the rough estimates of trend

D1t=(Yt-m1t)
j=1:33
wk=array(0)
for(k in 1:12){
wk[k]=mean(D1t[k+(j-1)*12])
}
length(wk)
sk=wk-mean(wk)  # This is the estimate of seasonal component

st=rep(sk,times=33)
plot(st,type="l")

D2t=Yt-st # Final deseasonalized data
plot(D2t,type="l",ylab = "dt",main="Final estimates of Trend Component from deseasonalised Data")

t=1:396
t1=t^2
m2t=predict(lm(D2t~t+t1))  # This is the estimate of Trend component
lines(t,m2t)


It=Yt-m2t-st
plot(It,type="l",ylab = "Ith", main = "Estimated Irregular Series")

## Testting for irregular component
n=396
u=0
for(i in 2:395)
{
  if(((It[i]>It[i+1])&&(It[i]>It[i-1]))||((It[i]<It[i+1])&&(It[i]<It[i-1])))
  {
    u=u+1
  }
}
u
e_u=2*(n-2)/3
v_u=(16*n-29)/90
z=(u-e_u)/v_u
z
qnorm(1-0.025)

# So it is accepted that It is a stationary process


### Checking whether It is WN or not
length(It)
sqrt(396)*cor(It[1:395],It[2:396])
qnorm(1-0.025)
# Reject H0 and conclude that It series is not WN

acf(It,main="Sample ACF plot of Estimated Irregular series")
pacf(It,main="Sample PACF plot of Estimated Irregular series")


###########
aic=array(0,dim=c(15,15))
for(i in 0:14)
{
  for(j in 0:14)
  {
    aic[(i+1),(j+1)]=0
    model=arima(It,order=c(i,0,j),method = "ML")
    aic[(i+1),(j+1)]=model$aic
  }
}
round(aic,2)
min(aic)

aic=matrix(aic,nrow=15,byrow=F)
rownames(aic)=c(0:14)
colnames(aic)=c(0:14)
aic

# Looking at aic matrix we choose p=13,q=12
p=13
q=12
model=arima(It,order=c(p,0,q),method = "ML")
model

### Forecasting Part

fc=predict(model,n.ahead=96,newxreg = NULL,se.fit = F)
fc  # predicted irregular values
sth = rep(sk,times = 8)
sth # predicted irregular values
mth = function(t){
  59.178120 + 0.223504*t - 0.000282*t^2
}
mth(397:492)
finfore=mth(397:492)+sth+fc
dd=seq(as.Date("2018-1-1"),as.Date("2025-12-1"),by="months")
length(dd)
aa=data.frame(dd,finfore)
colnames(aa)=c("Date", "  Forecasted Values")
View(aa)
write.csv(aa,"Forecasted.csv")
