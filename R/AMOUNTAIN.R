#' Illustration of weighted network simulation
#' 
#' Simulate a single weighted network
#' 
#' @param n number of nodes in the network
#' @param k number of nodes in the module, n < k
#' @param theta module node score follow the uniform distribution in range [theta,1]
#' 
#' @return a list containing network adjacency matrix, node score and module membership 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords simulation
#' 
#' @examples
#' pp <- networkSimulation(100,20,0.5)
#' moduleid <- pp[[3]]
#' netid <- 1:100
#' restp<- netid[-moduleid]
#' groupdesign=list(moduleid,restp)
#' names(groupdesign)=c('module','background')
#' \dontrun{library(qgraph)
#' pg<-qgraph(pp[[1]],groups=groupdesign,legend=TRUE)}
#' @export
#' 
networkSimulation<-function(n,k,theta){
    E <- matrix(runif(n*n, min=0, max=1),nrow = n)
    E[lower.tri(E)] = t(E)[lower.tri(E)]
    diag(E) = 1
    
    Z <- runif(n, min=0, max=1)
    moduleid <- sample(n, k)
    Z[moduleid] <- runif(k, min=theta, max=1)
    
    module <- matrix(runif(k*k, min=theta, max=1),nrow = k)
    module[lower.tri(module)] = t(module)[lower.tri(module)]
    diag(module) = 1
    
    E[moduleid,moduleid] <- module
    return (list(E,Z,moduleid))
}

#' Illustration of two-layer weighted network simulation
#' 
#' Simulate a two-layer weighted network
#' 
#' @param n1 number of nodes in the network1
#' @param k1 number of nodes in the module1, n1 < k1
#' @param theta1 module1 node score follow the uniform distribution in range [theta1,1]
#' @param n2 number of nodes in the network2
#' @param k2 number of nodes in the module2, n2 < k2
#' @param theta2 module2 node score follow the uniform distribution in range [theta2,1]
#' 
#' @return a list containing network1, network2 and a inter-layer links matrix 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @seealso \code{\link{networkSimulation}}
#' @keywords simulation
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
#' netid2 <- 1:n2
#' restp2<- netid2[-moduleid2]
#' ## labelling the groups
#' groupdesign=list(moduleid,restp,(moduleid2+n1),(restp2+n1))
#' names(groupdesign)=c('module1','background1','module2','background2')
#' twolayernet<-matrix(0,nrow=(n1+n2),ncol=(n1+n2))
#' twolayernet[1:n1,1:n1]<-pp[[1]]
#' twolayernet[(n1+1):(n1+n2),(n1+1):(n1+n2)]<-pp2[[1]]
#' twolayernet[1:n1,(n1+1):(n1+n2)] = A
#' twolayernet[(n1+1):(n1+n2),1:n1] = t(A)
#' \dontrun{library(qgraph)
#' g<-qgraph(twolayernet,groups=groupdesign,legend=TRUE)}
#' @export
#' 
twolayernetworkSimulation<-function(n1,k1,theta1,n2,k2,theta2){
    pp1<-networkSimulation(n1,k1,theta1)
    pp2<-networkSimulation(n2,k2,theta2)
    A <- matrix(0,nrow = n1,ncol = n2)
    E <-matrix(runif(k1*k2, min=0, max=1),nrow = k1)
    E[lower.tri(E)] = t(E)[lower.tri(E)]
    diag(E) = 0
    A[pp1[[3]],pp2[[3]]] <- E
    return (list(pp1,pp2,A))
}

#' Module Identification
#' 
#' Algorithm for Module Identification on single network
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
#' @seealso \code{\link{EuclideanProjectionENNORM}}
#' @keywords module identification
#' 
#' @examples
#' n = 100
#' k = 20
#' theta = 0.5
#' pp <- networkSimulation(n,k,theta)
#' moduleid <- pp[[3]]
#' ## use default parameters here
#' x <- moduleIdentificationGPFixSS(pp[[1]],pp[[2]],rep(1/n,n))
#' predictedid<-which(x[[2]]!=0)
#' recall <- length(intersect(predictedid,moduleid))/length(moduleid)
#' precise <- length(intersect(predictedid,moduleid))/length(predictedid)
#' Fscore <- (2*precise*recall/(precise+recall))
#' @export
#' 
moduleIdentificationGPFixSS <- function(W,z,x0,a=0.5,lambda=1,maxiter=1000){
    x = x0
    epsilon = 1e-6
    grad = -W%*%x-lambda*z
    f_x = -0.5*t(x)%*%W%*%x-lambda*(t(z)%*%x)
    func = numeric(maxiter)

    for (iteration in 1:maxiter){
            #y = x-1*grad
            #print(sum(y)+0.5*gamma*sum(y*y))
            func[iteration] = f_x
            #x_cand = EuclideanProjectionEN(x-1*grad,t=1,alpha = a)
            #x_cand = EuclideanProjection(x-1*grad,t=radius)
            x_cand = EuclideanProjectionENNORM (x-1*grad,t=1,alpha = a)
            if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
            x=x_cand
            grad = -W%*%x-lambda*z
            f_x = -0.5*t(x)%*%W%*%x-lambda*(t(z)%*%x)
            
    }
    return (list(func[1:(iteration-1)],x))
}

#' Euclidean projection on elastic net
#' 
#' Piecewise root finding algorithm for Euclidean projection on elastic net
#' 
#' @param y constant vector
#' @param t radius of elastic net ball
#' @param alpha parameter in elastic net: alpha x_1 + (1-alpha)*x_2^2=t
#' 
#' @return a list containing network adjacency matrix, node score and module membership 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references Gong, Pinghua, Kun Gai, and Changshui Zhang. "Efficient euclidean projections via piecewise root finding and its application in gradient projection." Neurocomputing 74.17 (2011): 2754-2766.
#' @keywords Euclidean projection
#' 
#' @examples
#' y=rnorm(100)
#' x=EuclideanProjectionENNORM(y,1,0.5)
#' sparistyx = sum(x==0)/100
#' @export
#' 
EuclideanProjectionENNORM <- function(y,t,alpha = 0.5){
    #f(theta) = \sum max(0,y_i-theta) - t
    maxiter = 100
    epsilon = 1e-9
    n = length(y)
    #theta_old = y[n]
    theta = 0
    x = numeric(n)
    if(alpha*sum(y)+(1-alpha)*sum(y*y) <= t){
        x = y
        return (x)
    } else {
        tmpx = (y-theta*alpha)/(1+(2-2*alpha)*theta)
        S = which(tmpx >= 0)
        a = (alpha^3-alpha^2)*length(S)-4*t*(1-alpha)^2
        b = -alpha^2*length(S)-4*t*(1-alpha)
        c = alpha*sum(y[S])+(1-alpha)*sum(y[S]*y[S])-t
        for (i in 1:maxiter) {
            theta = (-b-sqrt(b^2-4*a*c))/(2*a)
            # doesn't work for (-b+sqrt(b^2-4*a*c))/(2*a) why
            tmpx = (y-theta*alpha)/(1+(2-2*alpha)*theta)
            S = which(tmpx >= 0)
            a = (alpha^3-alpha^2)*length(S)-4*t*(1-alpha)^2
            b = -alpha^2*length(S)-4*t*(1-alpha)
            c = alpha*sum(y[S])+(1-alpha)*sum(y[S]*y[S])-t        
            if (a*theta^2+b*theta+c < epsilon)
                break
            #theta_old = theta
        }
    }
    for (i in 1:n) {
        x[i] = max(0,(y[i] - theta*alpha)/(1+(2-2*alpha)*theta))
    }
    return (x)
}

#' Module Identification for two-layer network
#' 
#' Algorithm for Module Identification on two-layer network
#' 
#' @param W1 edge score matrix of the network 1, n_1 x n_1 matrix
#' @param z1 node score vector of the network 1, n_1-length vector
#' @param x0 initial solution of network 1, n_1-length vector
#' @param W2 edge score matrix of the network 2, n_2 x n_2 matrix
#' @param z2 node score vector of the network 2, n_2-length vector
#' @param y0 initial solution of network 2, n_2-length vector
#' @param A inter-layer links weight, n_1 x n_2 matrix
#' @param lambda1 parameter in objective, coefficient of node score of network 1
#' @param lambda2 parameter in objective, coefficient of node score of network 2
#' @param lambda3 parameter in objective, coefficient of inter-layer links part
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing solution for network 1 and network 2 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{EuclideanProjectionENNORM}}
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
#' modres=moduleIdentificationGPFixSSTwolayer(pp[[1]],pp[[2]],rep(1/n1,n1),
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
moduleIdentificationGPFixSSTwolayer <- function(W1,z1,x0,W2,z2,y0,A,lambda1=1,
                                    lambda2=1,lambda3=1,maxiter=1000,a=0.5){
    x = x0
    y = y0
    epsilon = 1e-6
    grad_x = -W1%*%x-lambda1*z1-A%*%y
    grad_y = -W2%*%y-lambda2*z2-t(A)%*%x
    f_x = -0.5*t(x)%*%W1%*%x-lambda1*(t(z1)%*%x)-0.5*t(y)%*%W2%*%y-
        lambda2*(t(z2)%*%y)-t(x)%*%A%*%y
    func = numeric(maxiter)
    
    for (iteration in 1:maxiter){
        #y = x-1*grad
        #print(sum(y)+0.5*gamma*sum(y*y))
        func[iteration] = f_x
        x_cand = EuclideanProjectionENNORM(x-1*grad_x,t=1,alpha = a)
        #grad_x = -W1%*%x_cand-lambda1*z1-A%*%y
        
        y_cand = EuclideanProjectionENNORM(y-1*grad_y,t=1,alpha = a)
        #grad_y = -W2%*%y_cand-lambda2*z2-t(A)%*%grad_x
        #x_cand = EuclideanProjection(x-1*grad,t=radius)
        #x_cand = EuclideanProjectionENNORM (x-1*grad,t=1,alpha = a)
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon && 
               sum(abs(y_cand-y)^2)^(1/2) < epsilon){break}
        x = x_cand
        y = y_cand
        grad_x = -W1%*%x-lambda1*z1-A%*%y
        grad_y = -W2%*%y-lambda2*z2-t(A)%*%x
        
        f_x = -0.5*t(x)%*%W1%*%x-lambda1*(t(z1)%*%x)-0.5*t(y)%*%W2%*%y-
            lambda2*(t(z2)%*%y)-t(x)%*%A%*%y
    }
    return (list(x,y))
}

#' Module Identification for multi-layer network
#' 
#' Algorithm for Module Identification on multi-layer network sharing the same set of genes
#' 
#' @param W edge score matrix of the network, n x n matrix
#' @param zi node score vector of this layer, n-length vector
#' @param xrest consensus solution from the rest membership \code{\link{vecconsensus}}
#' @param x0 initial solution, n-length vector
#' @param lambdai parameter in objective, coefficient of node score of current layer
#' @param lambda parameter in objective, coefficient of node score of other layers
#' @param a parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}
#' @param maxiter maximal interation of whole procedure
#' 
#' @return a list containing solution for network 1 and network 2 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @references AMOUNTAIN
#' @seealso \code{\link{vecconsensus}}
#' @keywords module identification, two-layer
#' 
#' @examples
#' n=100
#' k=20
#' theta = 0.5
#' pp <- networkSimulation(n,k,theta)
#' W <- pp[[1]]
#' moduleid <- pp[[3]]
#' z1 <- runif(n, min=0, max=1)
#' z1[moduleid] <- runif(k, min=theta, max=1)
#' z2 <- runif(n, min=0, max=1)
#' z2[moduleid] <- runif(k, min=theta, max=1)
#' z3 <- runif(n, min=0, max=1)
#' z3[moduleid] <- runif(k, min=theta, max=1)
#' x0=rep(1/n,n)
#' x1=x0
#' x2=x0
#' x3=x0
#' xrest = vecconsensus(cbind(x2,x3))
#' x1 = moduleIdentificationGPFixSSMultilayer(W,z1,xrest,x0,lambdai=1,lambda=1,maxiter=100,a=0.5)
#' xrest = vecconsensus(cbind(x1,x3))
#' x2 = moduleIdentificationGPFixSSMultilayer(W,z2,xrest,x0,lambdai=1,lambda=1,maxiter=100,a=0.5)
#' xrest = vecconsensus(cbind(x1,x2))
#' x3 = moduleIdentificationGPFixSSMultilayer(W,z3,xrest,x0,lambdai=1,lambda=1,maxiter=100,a=0.5)
#' @export
#'
moduleIdentificationGPFixSSMultilayer <- function(W,zi,xrest,x0,lambdai=1,lambda=1,maxiter=1000,a=0.5){
    x = x0
    epsilon = 1e-6
    grad_x = -W%*%x-lambdai*zi-lambda*xrest
    
    f_x = -0.5*t(x)%*%W%*%x-lambdai*(t(zi)%*%x)-lambda*(t(x)%*%xrest)
    func = numeric(maxiter)
    
    for (iteration in 1:maxiter){
        func[iteration] = f_x[1,1]
        x_cand = EuclideanProjectionENNORM(x-1*grad_x,t=1,alpha = a)
        grad_x = -W%*%x-lambdai*zi-lambda*xrest
        if(sum(abs(x_cand-x)^2)^(1/2) < epsilon){break}
        x = x_cand
        grad_x = -W%*%x-lambdai*zi-lambda*xrest   
        f_x = -0.5*t(x)%*%W%*%x-lambdai*(t(zi)%*%x)-lambda*(t(x)%*%xrest)
    }
    return (x)
}

#' Vector consensus
#' 
#' a consensus of several vectors, keeping nonzero positions in all vectors
#' 
#' @param Mat membership matrix from other layers, n x l matrix
#' 
#' @return n-length vector 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' 
#' @export
#' 
vecconsensus <- function(Mat){
    n = dim(Mat)[1]
    p = c(length=n)
    for (i in 1:n){
        p[i]=prod(Mat[i,])
        if(p[i]!=0)
            p[i] = 1
    }
    return(p)
}