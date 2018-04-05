library(quadprog)

# my own version of the QP wrapper...needed to get rid of the 'stop' statements
my.solve.QP <- function (Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE) 
{
    n <- nrow(Dmat)
    q <- ncol(Amat)
    if (missing(bvec)) 
        bvec <- rep(0, q)
    if (n != ncol(Dmat)) 
        stop("Dmat is not symmetric!")
    if (n != length(dvec)) 
        stop("Dmat and dvec are incompatible!")
    if (n != nrow(Amat)) 
        stop("Amat and dvec are incompatible!")
    if (q != length(bvec)) 
        stop("Amat and bvec are incompatible!")
    if ((meq > q) || (meq < 0)) 
        stop("Value of meq is invalid!")
    iact <- rep(0, q)
    nact <- 0
    r <- min(n, q)
    sol <- rep(0, n)
    crval <- 0
    work <- rep(0, 2 * n + r * (r + 5)/2 + 2 * q + 1)
    iter <- rep(0, 2)
    res1 <- .Fortran("qpgen2", as.double(Dmat), dvec = as.double(dvec), 
        as.integer(n), as.integer(n), sol = as.double(sol), crval = as.double(crval), 
        as.double(Amat), as.double(bvec), as.integer(n), as.integer(q), 
        as.integer(meq), iact = as.integer(iact), nact = as.integer(nact), 
        iter = as.integer(iter), work = as.double(work), ierr = as.integer(factorized), 
        PACKAGE = "quadprog")
	#lhh modified the error handling so we can put the solver in a loop
    if (res1$ierr == 1) 
        list(status= -1, value = "constraints are inconsistent, no solution!")
    else if (res1$ierr == 2) 
        list(status=-2, value = "matrix D in quadratic function is not positive definite!")
	else
		list(status=0, solution = res1$sol, value = res1$crval, unconstrainted.solution = res1$dvec, 
			iterations = res1$iter, iact = res1$iact[1:res1$nact])
}


# convert the optimization problem into the standard form of the QP solver.
optimize <- function(df,b) {
	# We want to minimize (w,b): 1/2 w^T w
	# with the constraints:      X^T w >= b0
	# where xi^T w >= 1+b if yi = +1         (1)
	# where (-xi)^T w >= 1-b if yi = -1      (2)
	# NOTE; only the two dimensional case is implemented

	# number of rows in the dataframe
	n <- length(df$X1)

	# Construct the X matrix that holds the contraints
	# we need to multiply my the class label to make the constraints look the same
	# see equation (2) above
	x1 <- df$X1*df$Y
	x2 <- df$X2*df$Y
	Xmat <- matrix(c(x1,x2),2,n,byrow=TRUE)
	
	# Construct the b0 vector
	# we want the b0 vector to have 1-b entries for negative classes and 1+b for positive classes
	one <- rep(1,n) # construct a (1,1,1,...) vector
	b.temp <- rep(b,n) # construct a (b,b,b,....) vector
	b0 <- one + (df$Y*b) # construct a (1-b, 1+b, ....) vector

	# construct the remaining matrixes in order to make the problem fit into the 
	# standard QP framework
	Dmat       <- matrix(0,2,2)
	diag(Dmat) <- 1
	dvec       <- c(0,0)

	# call the solver
	my.solve.QP(Dmat,dvec,Xmat,bvec=b0)
}
	
# init the plot and place the dataset onto the plot
dataset <- function(df) {
	quartz(width=8,height=8)
	# setup the plot
	plot(0:8,0:8,type="n",main="Convex Optimization",xlab="X1",ylab="X2")
	#abline(h=0)
	#bline(v=0)
	
	# plot the classes
	for (i in 1:length(df$X1))
		if (df$Y[i] > 0 )
			points(df$X1[i],df$X2[i],col="red")
		else
			points(df$X1[i],df$X2[i],col="blue")
}

# plot the decision surface
surface <- function(w,b) {
	slope = -(w[1]/w[2])
	offset =  (b)/w[2]
	offset1 =  (b+1)/w[2]
	offset2 =  (b-1)/w[2]

	cat("slope = ", slope, "offset = ", offset,"\n")
	
	# plot the decision surface with supporting hyperplanes
	abline(offset,slope,lty="solid",lwd=2,col="green")
	abline(offset1,slope,lty="dashed")
	abline(offset2,slope,lty="dashed")
}

# compute the magnitude of a vector
size <- function(v) {
	sqrt(sum(v*v))
}


# we have to minimize over two parameters, the optimize takes care of the quadratic minimization
# here we do a simple grid search of an interval of b and use optimize to compute the w for each b.
search <- function(df,b.seq) {
	w.size.min <- 10000
	b.min <- 10000
	status <- -1

	for (b in b.seq) {
		soln <- optimize(df,b)
		if (soln$status == 0) {
			status = 0 # we found at least one solution
			w <- soln$solution
			w.size <- size(w)
			cat("b = ", b, "\t |w| = ", w.size,"\n")
			if (w.size < w.size.min) {
				w.size.min <- w.size
				w.min <- w
				b.min <- b
			}
		}
		else
			cat("b = ", b, "\t |w| = undefined\n")
	}
	
	if (status != 0)
		stop ("constraints inconsistent, no solution possible")
		
	list(w=w.min,b=b.min)
}

# driver program
run <- function() {
	# data set
	df <- as.data.frame(list("X1" = c(1,3,4,5,1,2,2.5,3), "X2"=c(6,7,1,3,4,1,1.5,1), "Y"=c(-1,-1,-1,-1,1,1,1,1)))
	# df <- as.data.frame(list("X1" = c(1,1,3,3,4), "X2"=c(2,4,4,1,2), "Y"=c(-1,-1,-1,1,1)))

	# start the plot
	dataset(df)

	# b value interval...[-20,20] with step .1
	b.window <- seq(-20,0,.1)
	# b.window <- seq(-20,20,.1)
	
	# search for the optimal values of w and b
	soln <- search(df,b.window)

	cat("w=",soln$w ,"b=", soln$b,"\n")
	surface(soln$w,soln$b)
}
