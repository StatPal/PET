library(PET)







































# abc <- test_vec(matrix(1:4, 2), 2, 3)
# print(abc)

# round(Gausfilter_test(3, 5), 3)
# round(Gausfilter_test(4, 5), 3)
# round(Gausfilter_test(5, 5), 3)
# round(Gausfilter_test(6, 5), 3)
# cat("\n")

# Gausfilter <- function(n, fwhm) {
#     twosparsq <- 2*((fwhm/2.3548)^2)
#     tmp2 <- exp(-((0:(n/2))^2)/twosparsq)  # tmp1 <- dnorm(0:(n/2), 0, spar)
#     tmp1 <- c(tmp2, rev(tmp2[1:(n/2)]))
#     tmp1/sum(tmp1)
# }

# Gausfilter_2 <- function(n, fwhm) {
#     spar <- fwhm/2.3548
#     tmp <- dnorm((-n/2+1):(n/2), 0, spar)
#     tmp/sum(tmp)
# }

# Gausfilter_2 <- function(n, fwhm) {
#     twosparsq <- 2*((fwhm/2.3548)^2)
#     tmp2 <- exp(-((0:(n/2))^2)/twosparsq)
#     tmp3 <- c(tmp2, rev(tmp2[1:(n/2)]))
#     tmp3/sum(tmp3)
# }

# Gausfilter_2 <- function(n, fwhm = 1, scaled= TRUE){
#     spar <- fwhm/2.3548
#     if (spar < 1e-16) {
#         ff <- rep(0,n)
#         ff[1] <- 1
#     } else{
#         ff <- exp(-((0:(n/2))/spar)^2/2)
#         ff <- c(ff, rev(ff[2: ((n + 1)/2)]))
#     }
#     if(scaled)
#         ff <- ff/sum(ff)
#     ff
# }


# round(Gausfilter_2(3, 5), 3)
# round(Gausfilter_2(4, 5), 3)
# round(Gausfilter_2(5, 5), 3)
# round(Gausfilter_2(6, 5), 3)

