## Now checking with the C code images.

library(PET)

P <- matrix(scan("PETcodes_Maitra/Hoff.128"), 128)
# P <- P[34:94,12:116]  ## Bug here??
# P <- P[32:100,12:116]
# dim(P)


rP_new <- newPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 128)
# rP_new <- markPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 128)$rData
rP <- matrix(scan("PETcodes_Maitra/Hoff.fwd"), 320, byrow=T)


cat('\n\n\n\n\n Backprojected images\n')
irP2_shrinked <- Backproj_R_shrinked(rP_new, nrow(P), ncol(P), mode="BF")
print(dim(irP2_shrinked$backfilter))
print(dim(irP2_shrinked$filter))


y <- matrix(scan("PETcodes_Maitra/y.dat"), 320)
Kty <- matrix(scan("PETcodes_Maitra/Kty.dat"), 256)  # padded size
eigens <- matrix(scan("PETcodes_Maitra/eigens_Maitra.dat"), 256)  # padded size


# irP2_orig <- BackProject_R_orig_dim(rP_new$rData, nrow(P), ncol(P), mode="BF")
# print(dim(irP2_orig$backfilter))
# print(dim(irP2_orig$filter))



# cat('\n\n\n\n\n\n')
irP3 <- iradon(rP_new, nrow(P), ncol(P), mode="BF")


PET:::SLSPRESS_2d(1.2, rP_new, nrow(P), ncol(P))
PET:::Optimize_SLS_2d(rP_new, nrow(P), ncol(P))


# final <- PET:::Smoothed_image(rP_new, nrow(P), ncol(P))





pdf('Rplots.pdf', width=12, height=6)
par(mfrow = c(1, 2))


cat("Check\n\n")


# Sinograms
rasterImage::rasterImage2(z=rP_new)
rasterImage::rasterImage2(z=rP)
# image(final)

range(rP_new)
range(rP)
# range(final)

## Backprojected
rasterImage::rasterImage2(z=irP2_shrinked$backfilter)
rasterImage::rasterImage2(z=Kty)

range(irP2_shrinked$backfilter)
range(Kty)

# filter
rasterImage::rasterImage2(z=irP2_shrinked$filter)
rasterImage::rasterImage2(z=1/sqrt(eigens))


range(irP2_shrinked$filter)
range(1/sqrt(eigens))
