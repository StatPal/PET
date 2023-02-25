library(PET)

P <- matrix(scan("PETcodes_Maitra/Hoff.128"), 128)
# P <- P[34:94,12:116]  ## Bug here??
P <- P[32:100,12:116]
dim(P)


# rP <- radon(P, mode="LI", ThetaSamples = 320, RhoSamples = 128)
rP_new <- markPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 135)

# par(mfrow = c(1, 2))
# image(rP$rData)
# image(rP_new$rData)





cat('\n\n\n\n\n\n')
irP2 <- Backproj_R(rP_new$rData, nrow(P), ncol(P), mode="BF")
irP2_shrinked <- Backproj_R_shrinked(rP_new$rData, nrow(P), ncol(P), mode="BF")
print(dim(irP2_shrinked$backfilter))
print(dim(irP2_shrinked$filter))

# irP2_orig <- BackProject_R_orig_dim(rP_new$rData, nrow(P), ncol(P), mode="BF")
# print(dim(irP2_orig$backfilter))
# print(dim(irP2_orig$filter))



# cat('\n\n\n\n\n\n')
irP3 <- iradon(rP_new$rData, nrow(P), ncol(P), mode="BF")


# PET:::SLSPRESS_2d(1.2, rP_new$rData, nrow(P), ncol(P))

# PET:::Optimize_SLS_2d(rP_new$rData, nrow(P), ncol(P))









pdf('Rplots.pdf', width=18, height=6)
par(mfrow = c(1, 3))

# Comparison between changed and the original algo.
image(irP2$irData)
image(irP2_shrinked$irData)
# image(irP2_orig$irData)
image(irP3$irData)
# Those are similar to iradon, so the code is working. 


# Visual comparison between image and backfilter
par(mfrow = c(1, 3))
image(irP2$irData)
image(irP2$backfilter)
image(irP2_shrinked$backfilter)
# image(irP2_orig$backfilter)
# image(irP2$filter)

dim(irP2$backfilter)
dim(irP2_shrinked$backfilter)


par(mfrow = c(1, 2))
dim(irP2$filter)
dim(irP2_shrinked$filter)
image(irP2$filter)
image(irP2_shrinked$filter)




# cat("\n\nAbs and Rel diff:\n")
# mean(abs(irP2$irData - irP3$irData))
# mean(abs(irP2$irData - irP3$irData))/mean(abs(irP2$irData))
# identical(irP2$irData, irP3$irData)


# dim(irP2$irData)


