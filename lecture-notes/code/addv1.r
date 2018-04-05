addv1 <- function(v) 
	{ 
		y <- c() 
		for (x in v) { 
			x1 <- x + 1 
			y <- c(y,x1)
		} 
		y
	}