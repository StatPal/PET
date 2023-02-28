## Chekcing for rectatngular matrixes,
# what happens when your original image is rectangular,
# Then the radon function and the markPossion does not match visually!
# Also, their definition of R and the Maitra's defn of N are diff. 
# Check that also. 

library(PET)

P <- matrix(scan("PETcodes_Maitra/Hoff.128"), 128)
## P <- P[34:94,12:116]  ## Bug here??
# P <- P[32:100,12:116]
# P <- P[20:100,12:116]
dim(P)


# rP <- radon(P, mode="LI", ThetaSamples = 320, RhoSamples = 128)
rP_new <- markPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 135)
rP <- newPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 135)

irP3 <- iradon(rP, nrow(P), ncol(P), mode="BF")
irP3_new <- iradon(rP_new$rData, nrow(P), ncol(P), mode="BF")

pdf('Rplots.pdf', width=18, height=6)
par(mfrow = c(1, 3))
image(P)
image(irP3$irData)
image(irP3_new$irData)





rP <- newPoisson(P, nSample=5e6, ThetaSamples = 320, RhoSamples = 135)

# cat('\n\n\n\n\n\n')
# irP2 <- Backproj_R(rP_new$rData, nrow(P), ncol(P), mode="BF")
# irP2_shrinked <- Backproj_R_shrinked(rP_new$rData, nrow(P), ncol(P), mode="BF")
# print(dim(irP2_shrinked$backfilter))
# print(dim(irP2_shrinked$filter))



# PET:::SLSPRESS_2d(1.2, rP_new$rData, nrow(P), ncol(P))
PET:::Optimize_SLS_2d(rP, nrow(P), ncol(P)) # nolint



