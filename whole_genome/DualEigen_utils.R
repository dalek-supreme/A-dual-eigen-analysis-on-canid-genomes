# principle.loadings <- function(vec,percentage=0.8,norm=norm.2.sqrt,criterion='norm') {
#     if (criterion == 'norm') {
#         sorted <- sort(vec)
#         v.norm <- norm(vec)
#         s <- array()
#         for (i in seq_along(vec)) {
#             s[i] <- vec[i]
#             if (norm(s)/v.norm>percentage)
#                 break 
#         }
#         return(order(vec)[seq(i)])
#     }
#     if (criterion == 'percentile') {
        
#     }
# }

# norm.2.sqare <- function(vec) {
#     n <- 0
#     for (i in seq_along(vec)){
#         n <- n + (vec[i])^2
#     }
#     return(n)
# }

# norm.2.sqrt <- function(v) { return(sqrt(norm.2.sqare(v)))}

# norm.p <- function(vec,p) {
#     n <- 0
#     for (i in seq_along(vec)){
#         n <- n + (vec[i])^p
#     }
#     return(abs(n)^(1/p))
# }

# norm.3 <- function(v) {return(norm.p(v,3))}

# quantile(vec,0.8)

principle.loadings <- function(vec,percentile=0.9,number,UseNumber=F) {
    if (!UseNumber){
        result <- vector()
        sval <- quantile(vec,1-percentile)
        lval <- quantile(vec,percentile)
        for (val in vec){
            if (val<sval | val>lval) {
                result<-append(result,val)
            }
        }
        return(sort(result))
    }
    else {
        sorted <- sort(vec)
        result <- vector()
        for (i in seq(min(number/2,length(vec)))) {
            result <- append(result,sorted[i])
        }
        for (i in length(vec)-seq(min(number/2,length(vec)))) {
            result <- append(result,sorted[i])
        }
        return(sort(result))
    }
}

non.principle.loadings <- function(vec,percentile=0.9,number,UseNumber=F) {
    if (!UseNumber){
        result <- vector()
        sval <- quantile(vec,1-percentile)
        lval <- quantile(vec,percentile)
        for (val in vec){
            if (val>sval & val<lval) {
                result<-append(result,val)
            }
        }
        return(sort(result))
    }
    else {
        sorted <- sort(vec)
        result <- sorted[seq(number/2,length(vec)-number/2)]
        return(result)
    }
}