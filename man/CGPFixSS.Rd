\name{CGPFixSS}
\alias{CGPFixSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Module Identification
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CGPFixSS(W, z, x0, a = 0.5, lambda = 1, maxiter = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{W}{edge score matrix of the network, n x n matrix}

\item{z}{node score vector of the network, n-length vector}

\item{x0}{initial solution, n-length vector}

\item{a}{parameter in elastic net the same as in \code{\link{EuclideanProjectionENNORM}}}

\item{lambda}{parameter in objective, coefficient of node score part}

\item{maxiter}{maximal interation of whole procedure}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a list containing function objective vector and the solution
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
AMOUNTAIN
}
\author{
%%  ~~who you are~~
Dong Li, \email{dxl466@cs.bham.ac.uk}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
n = 100
k = 20
theta = 0.5
pp <- networkSimulation(n,k,theta)
moduleid <- pp[[3]]
## use default parameters here
x <- CGPFixSS(pp[[1]],pp[[2]],rep(1/n,n))
predictedid<-which(x[[2]]!=0)
recall <- length(intersect(predictedid,moduleid))/length(moduleid)
precise <- length(intersect(predictedid,moduleid))/length(predictedid)
Fscore <- (2*precise*recall/(precise+recall))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ identification }% use one of  RShowDoc("KEYWORDS")
\keyword{ module }% __ONLY ONE__ keyword per line
