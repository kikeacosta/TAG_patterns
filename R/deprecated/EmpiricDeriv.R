emp.der <- function(x, y){
# simple function for "empirical" derivatives
    m <- length(x)
    a <- c(diff(y), 1)
    b <- c(diff(x), 1)
    ab <- a/b
    wei <- c(rep(1, m-1), 0)

    a1 <- c(1, -1*diff(y))
    b1 <- c(1, -1*diff(x))
    ab1 <- a1/b1
    wei1 <- c(0, rep(1, m-1))

    y1emp <- (ab*wei + ab1*wei1)/(wei+wei1)
    return(y1emp)
}

## for equally-spaced y
# x <- 1:10
# y <- runif(10)#c(2,3,4,5,6,2,3,4,5,6)
# y1emp <- emp.der(x,y)
# 
# A1 <- diff(diag(10),diff=1)
# A1 <- rbind(A1,0)
# A2 <- diff(diag(10),diff=1)
# A2 <- rbind(0, A2)
# 
# w <- rep(2,10)
# w[c(1,10)] <- 1
# A <- (A1+A2)/w
# y1 <- A%*%y
# 
# 
# cbind(y1, y1emp)
