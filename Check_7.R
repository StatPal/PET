## Checking Refft

## Now checking with the C code images.

library(PET)

P <- matrix(0, 64, 64)
P <- matrix(0, 32, 32)
P[1,1] <- 1
P[1,2] <- 2.5

n <- nrow(P)


orig <- matrix(scan("PETcodes_Maitra/refft.dat"), n, byrow=T)

image(PET:::ReFFTf2(P))
image(orig)


