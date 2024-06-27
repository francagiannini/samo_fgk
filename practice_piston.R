

#install.packages("sensobol")

library(sensobol)
library(tidyverse)

cycle_func <- function(i) {
  
  i <- sample_arrange[i,]
  
  A <- (i$P*i$S)+(19.62*i$M-(i$k*i$V0/i$S))
  
  V <- (i$S/2*i$k)*(sqrt((A^2)+(4*i$k*(i$P*i$V0/i$T0)*i$Ta))-A)
  
  C <- 2*pi*sqrt(i$M/(i$k+(i$S^2)*(i$P*V/i$T0)*i$Ta/(V^2)))
                 
  C
}

sample_arrange <-  sobol_matrices(
  matrices = c("A", "B", "AB"),
     N=1000,
     params = c("M","S","V0","k","P","Ta","T0")
   ) |> as_tibble()

sample_arrange$M <- qunif(sample_arrange$M, min = 30, max = 60)
sample_arrange$S <- qunif(sample_arrange$S, min = 0.005,max= 0.020)
sample_arrange$V0 <- qunif(sample_arrange$V0, min = 0.002, max = 0.010)
sample_arrange$k <- qunif(sample_arrange$k, min = 1000, max = 5000)
sample_arrange$P <- qunif(sample_arrange$P, min = 90000, max = 110000)
sample_arrange$Ta <- qunif(sample_arrange$Ta, min = 290, max = 296)
sample_arrange$T0 <- qunif(sample_arrange$T0, min = 340, max = 360)


output <- lapply(seq(1,nrow(sample_arrange),1), cycle_func)

output <- do.call(rbind,output)

hist(output)

sd(output)^2

index <-
  sobol_indices(
    Y = output,
    matrices = c("A", "B", "AB"),
    N = 1000,
    params = c("M", "S", "V0", "k", "P", "Ta", "T0")
  )

plot(index)

# Take A or B and compute the PCE simL@B

rm(list=ls())
library("randtoolbox")
source(paste0(getwd(),"/Benchmarks/piston.R"))
source(paste0(getwd(),"/BSPCE/BuildBSPCE_MAP.R"))

require(ggplot2)
Nvar <- 7
PCEType <- vector()
for(p in seq(1:Nvar)){
  PCEType <- c(PCEType,"Uniform")
}
PDF<-list()
PDF$Coeff <- matrix()
PDF$Coeff <- matrix(c(30 , 0.0125, 0.006 , 6.9, 9.0e+4 ,290, 296,
                      60 , 0.0025, 1/1500, 8.5, 10.0e+4,340, 360,
                      0.0, 0.0   , 0.0   ,0.0 , 11.0e+4,0.0, 0.0),
                    ncol=3,
                    nrow=Nvar)
PDF$Type <- vector()
PDF$Type <- c('Uniform','Normal','Laplace','LogUniform','Triangular','Uniform','Uniform')
## Randomized Sobol sequence (with digital shift)
#set.seed(271]
Nsample <- 2^7
U <- sobol(Nsample,Nvar)

X <- matrix(0,Nsample,Nvar)
#Transformation
X[,1] <- U[,1]*(PDF$Coeff[1,2]-PDF$Coeff[1,1])+PDF$Coeff[1,1]#Uniform
#
X[,2] <- qnorm(U[,2])*PDF$Coeff[2,2]+PDF$Coeff[2,1]#Normal
#
Ind1 <- which(U[,3]<0.5)
Ind2 <- which(U[,3]>=0.5)
X[Ind1,3] <- PDF$Coeff[3,1]+PDF$Coeff[3,2]*log(2*U[Ind1,3])
X[Ind2,3] <- PDF$Coeff[3,1]-PDF$Coeff[3,2]*log(2*(1-U[Ind2,3]))#Laplace
#
X[,4] <- exp(U[,4]*(PDF$Coeff[4,2]-PDF$Coeff[4,1])+PDF$Coeff[4,1])#LogUniform
#
Ind1 <- which(U[,5]<(PDF$Coeff[5,2]-PDF$Coeff[5,1])/(PDF$Coeff[5,3]-PDF$Coeff[5,1]))
Ind2 <- which(U[,5]>=(PDF$Coeff[5,2]-PDF$Coeff[5,1])/(PDF$Coeff[5,3]-PDF$Coeff[5,1]))
X[Ind1,5] <- PDF$Coeff[5,1]+sqrt(U[Ind1,5]*(PDF$Coeff[5,2]-PDF$Coeff[5,1])*(PDF$Coeff[5,3]-PDF$Coeff[5,1]))
X[Ind2,5] <- PDF$Coeff[5,3]-sqrt((PDF$Coeff[5,3]-PDF$Coeff[5,2])*(1-U[Ind2,5])*(PDF$Coeff[5,3]-PDF$Coeff[5,1]))#Triangle
#
X[,6] <- U[,6]*(PDF$Coeff[6,2]-PDF$Coeff[6,1])+PDF$Coeff[6,1]#Uniform
#
X[,7] <- U[,7]*(PDF$Coeff[7,2]-PDF$Coeff[7,1])+PDF$Coeff[7,1]#Uniform

#Model Call
y<-vector()
for(k in seq(from=1,to=Nsample)){
  y[k] <- piston(X[k,])
}
#PCE Building
PCE <- Build_SPCE(X,y)#
#Sensitivity indices Estimate
SA <- Compute_SI(PCE,X)#
AllIndices      <- SA
AllIndices[1:6] <- NULL

Epsilon=0.00
Abscisse <- seq(1:Nvar)
plot(Abscisse,SA$Si[,2], col = "red", pch=20, 
     xlab = "Input Number", ylab="Sobol' Indices",
     xlim=c(1,Nvar),ylim=c(min(SA$Si[,1]),(1+0.2)*max(SA$STi[,3])),xaxt="n")
axis(1,c(1:Nvar))
segments(Abscisse,SA$Si[,1],Abscisse,SA$Si[,3], col = "red")
segments(Abscisse-Epsilon,SA$Si[,1],Abscisse+Epsilon,SA$Si[,1], col = "red")
segments(Abscisse-Epsilon,SA$Si[,3],Abscisse+Epsilon,SA$Si[,3], col = "red")
points(Abscisse,SA$STi[,2], col = "blue", pch=21)
segments(Abscisse,SA$STi[,1],Abscisse,SA$STi[,3])
segments(Abscisse-Epsilon,SA$STi[,1],Abscisse+Epsilon,SA$STi[,1])
segments(Abscisse-Epsilon,SA$STi[,3],Abscisse+Epsilon,SA$STi[,3])
legend("topright",c("First-Order","Total-Order"), col=c("red","blue"),pch=c(20,21))

print("Unexplained amount of variance")
print(PCE$Res)
#The overall Indices
print("The overall estimated Sobol' indices")
print(unlist(AllIndices))
