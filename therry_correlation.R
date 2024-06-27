# practical Therry Thursday
d=3
N=1000
u1=runif(min=0,max=1,n=N)
u2=runif(min=0,max=1,n=N)
u3=runif(min=0,max=1,n=N)

u=cbind(u1,u2,u3)

#C+C´-eye(d)
C=matrix(c(1,-0.5,0,
         -0.5,1,0.8,
           0,0.8,1), nrow = 3)

Up=chol(C)

#sqrt(2)*erfinv(2*u(u_i,p)-1)
z=apply(u, 2, scale)

cov(z)
mean(z)

mu=seq(1,3,1) 

x=mu+(z%*%Up)

colnames(x) <- c("x1","x2","x3")

cov(x)

y=x[,"x1"]+x[,"x2"]+x[,"x3"]

# Go and check in simL@b and evaluate the differences

cbind(x,y)

cbind(u,y)
