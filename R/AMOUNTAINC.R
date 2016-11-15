# dyn.load(paste("../src/AMOUNTAIN", .Platform$dynlib.ext, sep = ""))

#' Module Identification
#' 
#' Call C version of \code{\link{moduleIdentificationGPFixSS }}
#' 
#' @param W edge score matrix of the network, n x n matrix
#' @param z node score vector of the network, n-length vector
#' @param x0 initial solution, n-length vector
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param lambda parameter in objective, coefficient of node score part
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing function objective vector and the solution 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{moduleIdentificationGPFixSS}}
#' @keywords module identification
#' 
#' @examples
#' n = 100
#' k = 20
#' theta = 0.5
#' pp <- networkSimulation(n,k,theta)
#' moduleid <- pp[[3]]
#' ## use default parameters here
#' x <- CGPFixSS(pp[[1]],pp[[2]],rep(1/n,n))
#' predictedid<-which(x[[2]]!=0)
#' recall <- length(intersect(predictedid,moduleid))/length(moduleid)
#' precise <- length(intersect(predictedid,moduleid))/length(predictedid)
#' Fscore <- (2*precise*recall/(precise+recall))
#' @export
#' 

CGPFixSS <- function(W,z,x0,a=0.5,lambda=1,maxiter=50){
	N = dim(W)[1]
	out <- .C("miGPFixSS",
    W = as.vector(W),
    z = as.vector(z),
    x0 = as.vector(x0),
    m_n = as.integer(N),
    x = as.vector(rep(1/N,N)),
    func = as.vector(rep(0,maxiter)),
    m_a = as.double(a),
    m_lambda = as.double(lambda),
    m_maxiter = as.integer(maxiter),
    m_protype = as.integer(1))
    return(list(out$func[1:(out$m_maxiter-1)],out$x))
}

#' Module Identification for two-layer network
#' 
#' Call C version of \code{\link{moduleIdentificationGPFixSSTwolayer}}
#' 
#' @param W1 edge score matrix of the network 1, n_1 x n_1 matrix
#' @param z1 node score vector of the network 1, n_1-length vector
#' @param x0 initial solution of network 1, n_1-length vector
#' @param W2 edge score matrix of the network 2, n_2 x n_2 matrix
#' @param z2 node score vector of the network 2, n_2-length vector
#' @param y0 initial solution of network 2, n_2-length vector
#' @param interlayerA inter-layer links weight, n_1 x n_2 matrix
#' @param lambda1 parameter in objective, coefficient of node score of network 1
#' @param lambda2 parameter in objective, coefficient of node score of network 2
#' @param lambda3 parameter in objective, coefficient of inter-layer links part
#' @param a1 parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param a2 parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing solution for network 1 and network 2 and objective
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{moduleIdentificationGPFixSSTwolayer}}
#' @keywords module identification, two-layer
#' 
#' @examples
#' n1=100
#' k1=20
#' theta1 = 0.5
#' n2=80
#' k2=10
#' theta2 = 0.5
#' ppresult <- twolayernetworkSimulation(n1,k1,theta1,n2,k2,theta2)
#' A <- ppresult[[3]]
#' pp <- ppresult[[1]]
#' moduleid <- pp[[3]]
#' netid <- 1:n1
#' restp<- netid[-moduleid]
#' pp2 <- ppresult[[2]]
#' moduleid2 <- pp2[[3]]
#' ## use default parameters here
#' modres=CGPFixSSTwolayer(pp[[1]],pp[[2]],rep(1/n1,n1),
#' pp2[[1]],pp2[[2]],rep(1/n2,n2),A)
#' predictedid<-which(modres[[1]]!=0)
#' recall = length(intersect(predictedid,moduleid))/length(moduleid)
#' precise = length(intersect(predictedid,moduleid))/length(predictedid)
#' F1 = 2*precise*recall/(precise+recall)
#' predictedid2<-which(modres[[2]]!=0)
#' recall2 = length(intersect(predictedid2,moduleid2))/length(moduleid2)
#' precise2 = length(intersect(predictedid2,moduleid2))/length(predictedid2)
#' F2 = 2*precise2*recall2/(precise2+recall2)
#' @export
#' 
CGPFixSSTwolayer <- function(W1,z1,x0,W2,z2,y0,interlayerA,
        lambda1=1,lambda2=1,lambda3=1,maxiter=100,a1=0.5,a2=0.5){
        N1 = dim(W1)[1]
        N2 = dim(W2)[1]
        out <- .C("miGPFixSSTwolayer",
	    W1 = as.vector(W1),
	    z1 = as.vector(z1),
	    x0 = as.vector(x0),
	    m_n1 = as.integer(N1),
	    x = as.vector(x0),
	    W2 = as.vector(W2),
	    z2 = as.vector(z2),
	    y0 = as.vector(rep(1/N2,N2)),
	    m_n2 = as.integer(N2),
	    y = as.vector(rep(1/N2,N2)),
	    A = as.vector(interlayerA),
	    func = as.vector(rep(0,maxiter)),
	    m_a1 = as.double(a1),
	    m_a2 = as.double(a2),
	    m_lambda1 = as.double(lambda1),
	    m_lambda2 = as.double(lambda2),
	    m_lambda3 = as.double(lambda3),
	    m_maxiter = as.integer(maxiter))
    return(list(out$x,out$y,out$func[1:out$m_maxiter]))
}


#' Module Identification for multi-layer network
#' 
#' Call C version of \code{\link{moduleIdentificationGPFixSSMultilayer}}
#' 
#' @param W edge score matrix of the network, n x n matrix
#' @param listzs a list of node score vectors, each layer has a n-length vector
#' @param x0 initial solution, n-length vector
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param lambda parameter in objective, coefficient of node score of other layers
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing solution for network 1 and network 2 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{moduleIdentificationGPFixSSMultilayer}}
#' @keywords module identification, multi-layer
#' 
#' @examples
CGPFixSSMultiLayer <- function(W,listzs,x0,a=0.5,lambda=1,maxiter=50){
	N = dim(W)[1]
	out <- .C("miGPFixSSMultilayer",
    W = as.vector(W),
    listz = as.vector(unlist(listzs)),
    m_L = as.integer(length(listzs)),
    x0 = as.vector(x0),
    m_n = as.integer(N),
    x = as.vector(rep(1/N,N)),
    func = as.vector(rep(0,maxiter)),
    m_a = as.double(a),
    m_lambda = as.double(lambda),
    m_maxiter = as.integer(maxiter))
    return(list(out$func[1:out$m_maxiter],out$x))
}
