rm(list = ls())
library(Matrix)
A=read.csv("mydata.csv")  # given data file is saved as mydata.csv
colnames(A)=c("x","y")
A
n=60
p=3
plot(A, main="Scatter plot of x(t) vs y(t)")
attach(A)
x
y


Q=function(beta){
  alpha2 = cov(exp(beta*x),y)/var(exp(beta*x))
  alpha1 = mean(y)-alpha2*mean(exp(beta*x))
  sum((y-alpha1-alpha2*exp(beta*x))^2)
}


beta=seq(-2,4, length=1000)
k=1
q=array(0)
for (i in beta) {
  q[k]=Q(i)
  k=k+1
}

plot(beta, q, main = "Plot the residual sum of squares as 
            a function of beta", ylab = "Qstar(beta)",type="l")


beta = beta[which.min(q)]
beta
alpha2 = cov(exp(beta*x),y)/var(exp(beta*x))
alpha2
alpha1 = mean(y)-alpha2*mean(exp(beta*x))
alpha1
abline(v=beta)


l2norm=function(v){
  sqrt(sum(v^2))
}
theta0=c(alpha1,alpha2,beta) # original initial choice

#  theta0=c(1,0.5,0.6) #testing initial choice
theta=c(0,0,0)
F.= matrix(0,nrow=n,ncol=p)
F.[,1]=1

iter=0
while (l2norm(theta-theta0)>1e-6 ) {

F.[,2]=exp(theta0[3]*x)
F.[,3]=theta0[2]*x*exp(theta0[3]*x)
r=y-theta0[2]*exp(theta0[3]*x)-theta0[1]
delta= solve(t(F.)%*%F.)%*%t(F.)%*%matrix(r)

theta=theta0
theta0=theta0+delta  

iter=iter+1
print(l2norm(theta-theta0))
}

theta
iter

v=theta[1]+theta[2]*exp(theta[3]*x)

plot(A, main="Graph showing observed and Fitted Responce")
lines(x,v,type="l")

F.[,2]=exp(theta[3]*x)
F.[,3]=theta[2]*x*exp(theta[3]*x)
F.
P=F.%*%solve(t(F.)%*%F.)%*%t(F.)
Y=matrix(y)
I=diag(c(rep(1,60)))
t(Y)%*%(I-P)%*%Y
t(F.)%*%F.
qf(1-0.05,3,n-3)
2.766438*0.05512061/19

v
sighs = sum((y-v)^2)/60 # sigma hat square
sighs

S11 = sighs*solve(t(F.)%*%F.)
S11
se = sqrt(diag(S11))

theta
se
CI=matrix(c(theta-qnorm(1-0.025)*se,theta+qnorm(1-0.025)*se),nrow = 3, byrow = F)
CI


2*sighs^2/n

r=y-v
plot(v,r,main="Residual Plot",xlab = "Fitted Responses", ylab = "Residuals")

shapiro.test(r)

qqnorm(r,main="QQPlot of Residuals")
qqline(r,col="red")

####Question 8
data=data.frame(x,y)
model1.nls=nls(y ~ a1 + a2*exp(b*x),data = data,
               start=list(a1 = 2.284847,a2 = 1.339692,b=1.489489),
               algorithm = "port")
summary(model1.nls)
s2=(0.0311^2)*57/60
s2
resi = residuals(model1.nls)
pred = predict(model1.nls)
plot(pred, resi, ylab = "Residuals", xlab = "Fitted Values", main = "Residual plot using nls")


### Question 9
theta0=c(2,1,0.1,1) 
theta=c(0,0,0,0)
F.= matrix(0,nrow=n,ncol=4)

iter=0
while (l2norm(theta-theta0)>1e-6 ) {
  
  F.[,1]=exp(theta0[3]*x)
  F.[,2]=exp(theta0[4]*x)
  F.[,3]=theta0[1]*x*exp(theta0[3]*x)
  F.[,4]=theta0[2]*x*exp(theta0[4]*x)
  
  r=y-theta0[1]*exp(theta0[3]*x)-theta0[2]*exp(theta0[4]*x)
  delta= solve(t(F.)%*%F.)%*%t(F.)%*%matrix(r)
  
  theta=theta0
  theta0=theta0+delta  
  
  iter=iter+1
  print(l2norm(theta-theta0))
}

theta
iter


data=data.frame(x,y)
model2.nls=nls(y ~ a1*exp(b1*x) + a2*exp(b2*x),data = data, 
            start=list(a1 = 1,a2 = 1,b1 = 0.1,b2 = 1))
summary(model2.nls)

model3.nls=nls(y ~ a1*exp(b1*x) + a2*exp(b2*x),data = data, 
               start=list(a1 = 1,a2 = 1,b1 = 0.1,b2 = 1),algorithm = "port" )
summary(model3.nls)

plot(A, main="Graph showing observed and Fitted Responce")
v=theta[1]*exp(theta[3]*x)+theta[2]*exp(theta[4]*x)
lines(x,v,type="l",col="red")

sum((y-v)^2)/60 # sigma hat square

r=y-v
plot(v,r,main="Residual Plot",xlab = "Fitted Responses", ylab = "Residuals")
shapiro.test(r)
qqnorm(r)
qqline(r,col="red")
