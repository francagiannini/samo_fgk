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

sobol