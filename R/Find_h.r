Gausfilter <- function(n, fwhm) {
    spar <- fwhm/2.3548
    if (spar < 1e-16) {
        tmp <- rep(0, n)
        tmp[1] <- 1
    } else {
        tmp <- exp(-((0:(n / 2)) / spar)^2 / 2)
        tmp <- c(tmp, rev(tmp[2: ((n + 1) / 2)]))
    }
    tmp/sum(tmp)
}
## Doing same cacluation again - change it later

# Only even cases
Gausfilter_2d <- function(m, n, fwhm) {
    g <- Gausfilter(m, fwhm)
    h <- Gausfilter(n, fwhm)
    return(g %*% t(h))
}

fftGauss_2d <- function(m, n, fwhm) {
    f <- Gausfilter_2d(m, n, fwhm)
    dum <- f[1, 1]
    f <- PET:::ReFFTf2(f)
    f * m * n * dum / sum(f)
}


convolve_2d <- function(f, g, pre.plan=NULL) {

    cf <- fft(f, inverse = FALSE)
    cg <- fft(g, inverse = FALSE)

    cf <- cf * cg * sqrt(nrow(cf) * nrow(cg))

    cf <- fft(cf, inverse = TRUE)
    cf <- cf / sqrt(length(cf))
    return(Re(cf))
}




# f -> K(V(D^{-1}f))
U1f_PET <- function(f, filter, ThetaSamples, RhoSamples){
    ff <- f * sqrt(filter)      # D^{-1}  # BUG, filter is multiplied, not divided
    ff <- PET:::ReFFTi2(ff)      # V

    # ThetaSamples <- nrow(f)
    # RhoSamples <- ncol(f)
    ## This was the bug. Those were the m and n values

    tmp <- PET:::radon(ff, mode="LI", ThetaSamples = ThetaSamples, RhoSamples = RhoSamples)$rData  # K
    return(tmp)
}



predictLS <- function(y, XSamples, YSamples)
# For any corrected data y, this function calculates the forward projection of
#    the LS reconstruction. The returned value is y'Ff where Ff is the projected
#    LS reconstruction.
{
    ThetaSamples <- nrow(y)
    RhoSamples <- ncol(y)

    f <- PET:::Backproj_R(y, XSamples, YSamples)
    Ff <- PET:::radon(f$irData, mode="LI", ThetaSamples = ThetaSamples, RhoSamples = RhoSamples)$rData

    sum(y*Ff)
}



# Basically modifies the Z1 and outputs the Z2 value.
# y is the sinogram
Uty_PET <- function(y, XSamples, YSamples){  # 320x135, 69, 105
    # z1 = U_1'y = D^{-1}V'K'y
    tmp <- PET:::Backproj_R_shrinked(y, XSamples, YSamples)

    # cat('\n\nAfter backproject\n')
    # print(dim(y))               # 320 x 135
    # print(XSamples)             # 69
    # print(dim(tmp$backfilter))  # 69 x 105  # correct
    # print(dim(tmp$filter))      # 69 x 105  # correct

    Kty <- tmp$backfilter           # K'
    Uty <- PET:::ReFFTf2(Kty)       # V'
    z1 <- Uty*sqrt(tmp$filter)      # D^{-1}, filter multiplied, not divided
    # print(any(tmp$filter<=0))

    # print('dim(z1)')
    # print(dim(z1))              # 69 x 105 # correct
    sum2 <- sum(z1^2)
    z2tz2 <- sum(y^2) - sum2


    # ## Next correction due to some reason
    U1f <- U1f_PET(z1, tmp$filter, nrow(y), ncol(y))  # Same dim as z1
    # cat('\n\nHaha\n')
    # print(dim(z1))  # 69 x 105
    # print(dim(y))   # 320 x 135
    # print(dim(U1f)) # 69 x 105  -- WRONG, corrected now
    sum3 <- sum(y*U1f)
    sum4 <- predictLS(y, XSamples, YSamples)
    dum <- sqrt(sum4/sum2)
    z1 <- z1 * dum;
    z2tz2 <- sum(y^2) - sum4

    return(list(z1 = z1, z2tz2 = z2tz2, filter = tmp$filter))
}

SLSPRESS_2d <- function(fwhm, y, XSamples, YSamples) {
    if (!(is.matrix(y)))
      stop("'rData' has to be of type 'matrix'.")


    if(fwhm < 0){
        return(Inf)
    } else {
        mma <- nrow(y)  # 320
        mmd <- ncol(y)  # 135
        tmp <- Uty_PET(y, XSamples, YSamples)

        m <- nrow(tmp$z1)  # CHECK
        n <- ncol(tmp$z1)

        eigfwhm <- fftGauss_2d(m, n, fwhm)
        trfwhm <- sum(eigfwhm) / (mma * mmd - XSamples * YSamples)
        sum_val <- sum( ( (1-eigfwhm) * tmp$z1)^2 )

        # Why that is named squared in the original file???
        sum_val <- sum_val + tmp$z2tz2 * (1+2.*trfwhm+trfwhm*trfwhm)
        return(sum_val)
    }
}

Optimize_SLS_2d <- function(y, XSamples, YSamples){

    if (!(is.matrix(y)))
      stop("'rData' has to be of type 'matrix'.")
    abstol = 1e-12; lower = 1e-6; upper = 20;

    # cat("\nExample values\n")
    # print(SLSPRESS_2d(0, y=y, XSamples=XSamples, YSamples=YSamples))
    # print(SLSPRESS_2d(1, y=y, XSamples=XSamples, YSamples=YSamples))
    # print(SLSPRESS_2d(10, y=y, XSamples=XSamples, YSamples=YSamples))
    # print(SLSPRESS_2d(20, y=y, XSamples=XSamples, YSamples=YSamples))
    ## Possible bug check

    fnl <- optim(1.0, SLSPRESS_2d, method="Brent", 
                lower=lower, upper=upper, 
                y=y, XSamples=XSamples, YSamples=YSamples)
    cat("\n\nmin = %f at %f\n", fnl$val, fnl$par)
    return(fnl$par)
}


Smoothed_image <- function(y, XSamples, YSamples){
    tmp <- PET:::Backproj_R_shrinked(y, XSamples, YSamples)
    fwhm <- Optimize_SLS_2d(y, XSamples, YSamples)
    print(fwhm)
    mat <- Gausfilter_2d(XSamples, YSamples, fwhm)

    return(convolve_2d(mat, tmp$irData))
}



# Step_3a <- function(K, lambda_tilde, filter){
    
    
#     ## Use FFT to get V' lambda_tilde 
#     V_lambda = 

#     ## z_1
#     z1 = solve(D_dot, V_lambda)
# }

# Step_3b <- function(y, z1){
#     z2z2 = sum(y * y) - sum(z1 * z1)

#     Omega_h = diag(S_h_eig(h))

#     zeta = t(z1) * (diag() - Omega_h)^2 * z1 + # optimize
#            z2z2 * sum(Omega_h)/(n-p)
    
#     return(zeta)
# }

