test_vec <- function(a, b, val){
    print(5)
    ab <- .C("test_vec_C",
       a = as.double(a),
       b = as.double(b), 
       as.integer(length(a)), 
       as.double(val), 
       PACKAGE="PET")
    # print(ab$a)
    # print(ab$b)
    return(list(ab$a, ab$b))
}


Gausfilter_test <- function(n, fhwm){
    ab <- .C("Gausfilter_test",
       as.integer(n),
       as.double(fhwm), 
       f = as.double(rep(0, n)), 
       PACKAGE="PET")$f
    return(ab)
}
