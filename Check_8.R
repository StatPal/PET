## Unit case test

## Now checking with the C code images.

library(PET)
set.seed(1)

P <- matrix(0, 128, 128)
# P[1,127] <- 1

P[128,128] <- 1



P <- matrix(0, 32, 32)
# P[16,16] <- 1
P[1,1] <- 1

n <- nrow(P)

# # Sq matrix
# unpad <- function(A){
# 	n <- nrow(A)/2
# 	A[(n/2):(3*n/2), (n/2):(3*n/2)]
# }

# orig <- matrix(scan("PETcodes_Maitra/orig.dat"), n)
# y <- matrix(scan("PETcodes_Maitra/y.dat"), 320, byrow = T)  ## Not 320
# # y <- matrix(scan("PETcodes_Maitra/Hoff.fwd"), 320, byrow=T)
# # y <- matrix(scan("PETcodes_Maitra/Hoff.fwdnonoise"), 128)
# Kty <- matrix(scan("PETcodes_Maitra/Kty.dat"), 2*n)  # padded size
# eigens <- matrix(scan("PETcodes_Maitra/eigens_Maitra.dat"), 2*n)  # padded size
# final <- matrix(scan("PETcodes_Maitra/Hoff.fwd"), n)
# final <- final/sum(final)




# rP_new_old <- newPoisson(P, nSample=1e7, ThetaSamples = 320, RhoSamples = 128, mode="LI")
# write(t(rP_new), file="PETcodes_Maitra/rP_new.dat")
# rP_new <- matrix(scan("PETcodes_Maitra/rP_new.dat"), 320, byrow = T)
# dim(rP_new)



# rP_new_old <- radon(P, ThetaSamples = 320, RhoSamples = 128, mode="LI")$rData
# # rP_new_old <- newPoisson(P, nSample=1e7, ThetaSamples = 320, RhoSamples = 128, mode="LI")

# rP_new_old <- PET:::Forward_R(P, ThetaSamples = 320, RhoSamples = 128, mode="LI")$rData




#rP_new_old <- radon(P, ThetaSamples = 10, RhoSamples = 5, mode="LI")$rData
# rP_new_old <- PET:::Forward_R(P, ThetaSamples = 10, RhoSamples = 5, mode="LI")$rData


rP_new <- radon(P, ThetaSamples = 320, RhoSamples = 128, mode="LI")$rData
rP_new_old <- PET:::Forward_R(P, ThetaSamples = 320, RhoSamples = 128, mode="LI", RhoMin=0, XYmin = 0, DeltaRho = 1/128)$rData
## Original code does not work



# # rP_new <- y
# rP_new <- rP_new_old
class(rP_new)

image(rP_new)
image(rP_new_old)

rP_new <- t(rP_new_old)


# stop()


# cat('\n\n\n\n\n Backprojected images\n')
# irP2_shrinked <- Backproj_R_shrinked(rP_new, nrow(P), ncol(P), mode="BF", 
# 					DeltaX=sqrt(.50)*(ncol(rP_new)/nrow(P)), 
#                    	DeltaY=sqrt(.50)*(ncol(rP_new)/ncol(P)))

# print(max(irP2_shrinked$filter))

# # cat('\n\n\n\n\n\n')
# irP3 <- iradon(rP_new, nrow(P), ncol(P), mode="BF")


# # irP4 <- PET:::Smoothed_image(rP_new, nrow(P), ncol(P))




# pdf('Rplots.pdf', width=12, height=6)
# par(mfrow = c(1, 2))


# cat("Check\n\n")

# # Original Image
# rasterImage::rasterImage2(z=P, main = "original image")
# rasterImage::rasterImage2(z=orig)
# # rasterImage::rasterImage2(z=P)

# # Sinograms
# rasterImage::rasterImage2(z=rP_new, main = "Sinograms")
# rasterImage::rasterImage2(z=y)

# ## Backprojected
# rasterImage::rasterImage2(z=irP2_shrinked$backfilter, main = "backprojected Image")
# rasterImage::rasterImage2(z=unpad(Kty))


# ## FFT of filter
# FFT_backproj <- irP2_shrinked$backfilter
# rasterImage::rasterImage2(z=PET:::ReFFTf2(FFT_backproj), main = "FFT of backproj")
# rasterImage::rasterImage2(z=PET:::ReFFTf2(Kty))
# rasterImage::rasterImage2(z=PET:::ReFFTi2(FFT_backproj), main = "iFFT of backproj")
# rasterImage::rasterImage2(z=PET:::ReFFTf2(Kty))


# # filter
# rasterImage::rasterImage2(z=irP2_shrinked$filter, main = "Filter (Fourier domain possibly)")
# rasterImage::rasterImage2(z=unpad(1/sqrt(eigens)))
# range(irP2_shrinked$filter)



# ## FFT of filter
# FFT_filter <- irP2_shrinked$filter
# FFT_filter_2 <- 1/sqrt(eigens)
# rasterImage::rasterImage2(z=PET:::ReFFTf2(FFT_filter), main = "FFT of Filter")
# rasterImage::rasterImage2(z=PET:::ReFFTf2(FFT_filter_2))
# rasterImage::rasterImage2(z=PET:::ReFFTi2(FFT_filter), main = "iFFT of Filter")
# rasterImage::rasterImage2(z=PET:::ReFFTi2(FFT_filter_2))


# abc <- PET:::ReFFTi2(irP2_shrinked$filter)
# # abc <- PET:::ReFFTi2(array(0.1, dim=dim(irP2_shrinked$filter)))
# abc_2 <- PET:::ReFFTi2(FFT_filter_2)

# rasterImage::rasterImage2(z=PET:::ReFFTf2(RFASTfMRI::Gauss.smooth(abc, c(.1, .1))), main = "iFFT of Filter + Smooth + FFT")
# rasterImage::rasterImage2(z=PET:::ReFFTf2(RFASTfMRI::Gauss.smooth(abc_2, c(.1, .1))))


# # Final
# rasterImage::rasterImage2(z=irP2_shrinked$irData/sum(irP2_shrinked$irData), main = "Final estimate (without smoothing)")
# # rasterImage::rasterImage2(z=irP3$irData)
# rasterImage::rasterImage2(z=final, main = "Final estimate from previous code without smoothing")
# # rasterImage::rasterImage2(z=irP4)


# which(irP2_shrinked$irData==max(irP2_shrinked$irData), arr.ind = T)
# which(final==max(final), arr.ind = T)

# irP2_shrinked$irData[80:84,80:84]



irP3 <- iradon(rP_new, nrow(P), ncol(P), mode="BF")
BF_method <- which(irP3$irData==max(irP3$irData), arr.ind = T)

irP3 <- iradon(rP_new, nrow(P), ncol(P), mode="FB")
FB_method <- which(irP3$irData==max(irP3$irData), arr.ind = T)

irP3 <- iradonIT(rP_new, nrow(P), ncol(P), mode="EM")
EM_method <- which(irP3$irData==max(irP3$irData), arr.ind = T)


rbind(BF_method, FB_method, EM_method)
