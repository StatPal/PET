FFTw <- function(x, pre.plan = NULL, inverse=FALSE)
{
    ## calculates the forward and backward Fast Fourier Transform of a complex 
    ## 1-dimensional sequence (using FFTW, the Fastest Fourier Transform of
    ## the West). 
    ## pre.plan is what speeds up FFTW, as I understand it.

    if (is.null(pre.plan)) 
        pre.plan <- fftw::planFFT(n = length(x))

    fc <- fftw::FFT(x, plan = pre.plan, inverse = inverse) 

    fc <- fc/sqrt(length(fc))    
    fc
}

ReFFTf <- function(x, pre.plan = NULL) {
    ## calculates the forward transform of a given vector x (that is not
    ## necessarily symmetric). Note that the input sequence is replaced by its
    ## transform on return.  We transform x -> W'x where W is a real matrix
    ## diagonalizing  a symmetric circulant matrix.

    b <- FFTw(x, pre.plan = pre.plan, inverse = FALSE)

    n <- length(x)
    f <- vector(length = n)
    f[1] <- Re(b[1])
    
    if ((n %% 2)==0)
        f[n/2 + 1]=Re(b[n/2+1]);    
    ## These are the (one/two for n odd/even) real eigenvector(s) of the  complex Fourier matrix anyway.

    f[2:ceiling(n/2)] <- Re(b[2:ceiling(n/2)] + b[n:floor(n/2+2)])/sqrt(2)
    ##second quadrant
    f[n:floor(n/2+2)] <- Im(b[2:ceiling(n/2)] - b[n:floor(n/2+2)])/sqrt(2)
    f
}


ReFFTi <- function(x, pre.plan = NULL) {
    ## calculates the backward transform of a given vector f (not necessarily
    ## symmetric). Note that the input sequence is replaced by its transform on
    ## return. We transform f -> Wf, where W is a real matrix diagonalizing a 
    ## symmetric circulant matrix. 

    n <- length(x)
    zcf  <- rep(0, n) -> zsf
    zcf[1:(floor(n/2) + 1)] <- x[1:(floor(n/2) + 1)]
    zsf[(floor(n/2) + 2):n] <- -x[(floor(n/2) + 2):n]

    zcf[1] <- zcf[1]/sqrt(2)
    if ((n %% 2)==0)  # n even
        zcf[n/2 + 1]  <- zcf[n/2 + 1]/sqrt(2); 

    zs <- FFTw(zsf, pre.plan = pre.plan)*sqrt(2) ## not inverse here 
    zc <- FFTw(zcf, pre.plan = pre.plan)*sqrt(2)

    f <- (Re(zc)+Im(zs))
    f
}

ReFFTf2 <- function(x, pre.plan.1 = NULL, pre.plan.2 = NULL) {
    ## x is 2D array, pre.plan.1 is the plan in the first dimension, and
    ## pre.plan.2 is the plan in the second dimension.

    n <- dim(x)[1]
    m <- dim(x)[2]

    if (is.null(pre.plan.1))
        pre.plan.1 <- fftw::planFFT(n = n)
    if (is.null(pre.plan.2))
        pre.plan.2 <- fftw::planFFT(n = m)

    x.refft <- apply(X = x, MARGIN = 1, FUN = ReFFTf, pre.plan = pre.plan.2)
    ## note that the dimension of x.refft is transpose that of x
    
    apply(X=x.refft, MARGIN = 1, FUN = ReFFTf, pre.plan = pre.plan.1)
}


ReFFTi2 <- function(x, pre.plan.1 = NULL, pre.plan.2 = NULL) {
    ## x is 2D array, pre.plan.1 is the plan in the first dimension, and
    ## pre.plan.2 is the plan in the second dimension.

    n <- dim(x)[1]
    m <- dim(x)[2]

    if (is.null(pre.plan.1))
        pre.plan.1 <- fftw::planFFT(n = n)
    if (is.null(pre.plan.2))
        pre.plan.2 <- fftw::planFFT(n = m)

    x.refft <- apply(X = x, MARGIN = 1, FUN = ReFFTi, pre.plan = pre.plan.2)
    ## note that the dimension of x.refft is transpose that of x
    apply(X=x.refft, MARGIN = 1, FUN = ReFFTi, pre.plan = pre.plan.1)
}


freq_gen <- function(ms) {
    rr <- lapply(ms,function(m){
        q <- floor(m/2)
        if(m %% 2 == 0){
            xx <- c(0:q,(q-1):1)/q
        } else{
            xx <- c(0:q,q:1)/q
        }
        xx
    })
    ones <- lapply(ms, function(m) rep(1,m))
    toret <- array(0,ms)
    for(k in 1:length(ms)){
        retd <- ones
        retd[[k]] <- rr[[k]]
        toret <- toret + Reduce("%o%",retd)^2
    }
    sqrt(toret)
}

## yes, so if you have a 4000 x 6000 matrix Y and you want only 1/10 of the smallest frequencies you do
## xx <- freq_gen(c(4000,6000))
## yy[xx<1/10]
## xx <- freq_gen(c(4000,6000))
## Usage: Y[xx<1/10] where Y is the output from ReFFT2
## but the cool thing is that if you have a 188,433- valued vector yy containing 1/10 of the frequencies, you can extract 1/20 of the frequencies by doing yy[ which(xx<1/10) %in% which(xx<1/20)]
