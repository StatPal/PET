## Unit case test

## Now checking with the C code images.

library(PET)

P <- matrix(0, 64, 64)
P <- matrix(0, 128, 128)
# P <- matrix(0, 32, 32)
P[1,1] <- 1

n <- nrow(P)

rP_new <- newPoisson(P, nSample=5e5, ThetaSamples = 320, RhoSamples = 128, mode="NN")

cat('\n\n\n\n\n Backprojected images\n')
irP2_shrinked <- Backproj_R_shrinked(rP_new, nrow(P), ncol(P), mode="BF")


# cat('\n\n\n\n\n\n')
irP3 <- iradon(rP_new, nrow(P), ncol(P), mode="BF")


irP4 <- PET:::Smoothed_image(rP_new, nrow(P), ncol(P))




pdf('Rplots.pdf', width=12, height=6)
par(mfrow = c(1, 2))

orig <- matrix(scan("PETcodes_Maitra/orig.dat"), n)
y <- matrix(scan("PETcodes_Maitra/y.dat"), 128)  ## Not 320
y <- matrix(scan("PETcodes_Maitra/Hoff.fwd"), 128)
# y <- matrix(scan("PETcodes_Maitra/Hoff.fwdnonoise"), 128)
Kty <- matrix(scan("PETcodes_Maitra/Kty.dat"), 2*n)  # padded size
eigens <- matrix(scan("PETcodes_Maitra/eigens_Maitra.dat"), 2*n)  # padded size
final <- matrix(scan("PETcodes_Maitra/Hoff.GCVP"), n)
final <- final/sum(final)


cat("Check\n\n")

# Original Image
rasterImage::rasterImage2(z=P, main = "original image")
rasterImage::rasterImage2(z=orig)
# rasterImage::rasterImage2(z=P)

# Sinograms
rasterImage::rasterImage2(z=rP_new, main = "Sinograms")
rasterImage::rasterImage2(z=y)

## Backprojected
rasterImage::rasterImage2(z=irP2_shrinked$backfilter, main = "backprojected Image")
rasterImage::rasterImage2(z=Kty)


## FFT of filter
FFT_backproj <- irP2_shrinked$backfilter
rasterImage::rasterImage2(z=PET:::ReFFTf2(FFT_backproj), main = "FFT of backproj")
rasterImage::rasterImage2(z=PET:::ReFFTf2(Kty))
rasterImage::rasterImage2(z=PET:::ReFFTi2(FFT_backproj), main = "iFFT of backproj")
rasterImage::rasterImage2(z=PET:::ReFFTf2(Kty))


# filter
rasterImage::rasterImage2(z=irP2_shrinked$filter, main = "Filter (Fourier domain possibly)")
rasterImage::rasterImage2(z=1/sqrt(eigens))
range(irP2_shrinked$filter)



## FFT of filter
FFT_filter <- irP2_shrinked$filter
rasterImage::rasterImage2(z=PET:::ReFFTf2(FFT_filter), main = "FFT of Filter")
rasterImage::rasterImage2(z=PET:::ReFFTf2(1/sqrt(eigens)))
rasterImage::rasterImage2(z=PET:::ReFFTi2(FFT_filter), main = "iFFT of Filter")
rasterImage::rasterImage2(z=PET:::ReFFTf2(1/sqrt(eigens)))


# Final
rasterImage::rasterImage2(z=irP2_shrinked$irData/sum(irP2_shrinked$irData), main = "Final estimate (without smoothing)")
# rasterImage::rasterImage2(z=irP3$irData)
rasterImage::rasterImage2(z=final, main = "Final estimate from previous code")
# rasterImage::rasterImage2(z=irP4)