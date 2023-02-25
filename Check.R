library(PET)
# abc <- test_vec(matrix(1:4, 2), 2, 3)
# print(abc)

# Gausfilter_test(5, 0.5)


# P <- phantom(n=128, design="B")
P <- matrix(scan("PETcodes_Maitra/Hoff.128"), 128)
dim(P)

pdf('Rplots.pdf', width=12, height=6)
par(mfrow = c(1, 2))
# rP_new <- radon(P, mode="LI", ThetaSamples = 320, RhoSamples = 128)
# image(rP_new$rData)
# dim(rP_new$rData)

# rPnew <- matrix(scan("PETcodes_Maitra/Hoff.fwdnonoise"), ncol=128, nrow=320, byrow=T)
# image(rPnew)
# dim(rPnew)





cat("\n\nNow with noise \n")
# Checking sinograms
P1 <- P[12:116,12:116]
rP <- markPoisson(P1, nSample=5e5, ThetaSamples = 320, RhoSamples = 128)
print(dim(rP$rData))
# image(rP$rData)

# rPnew <- matrix(scan("PETcodes_Maitra/Hoff.fwd"), ncol=128, nrow=320, byrow=T)
# dim(rPnew)
# image(rPnew)

# # rPnew <- rP$rData

irP2 <- Backproj_R(rP$rData, nrow(P1), ncol(P1), mode="BF")


# # irP1 <- iradon(rPnew, nrow(P), ncol(P), mode="FB", FilterTyp = 'Ramp') # Hamming filter was using inbuilt smoothing.
irP2 <- iradon(rP$rData, nrow(P1), ncol(P1), mode="BF")
# image(P)
# image(irP2$irData)


A <- matrix(scan('InvMyImage_after_backproject.dat'), 128, 128, byrow=T)
dim(P1)
dim(A)
dim(irP2$irData)

image(P1)
image(A)

# cat("\n\n\n\n\n\n\n")
# irP3 <- iradon_smoothed(rPnew, nrow(P), ncol(P))

# # pdf('Rplots.pdf', width = 12, height=3)
# # par(mfrow=c(1,4))
# # image(P)
# # image(irP1$irData)
# image(irP2$irData)
# # image(irP3$irData)

# # remove.packages("TestPack", lib="~/R/x86_64-pc-linux-gnu-library/4.2")
# # remove.packages("PET", lib="~/R/x86_64-pc-linux-gnu-library/4.2")



# ## First check the forward transformation as that seems different in breadth


