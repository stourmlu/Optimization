source('Optimization.R');

########## PROBLEM: minimize banana function, starting from x0 as defined below ##########
x0 <- c(-1.9, 2);

banana_func <- function(x, returnGrad=FALSE, returnHess=FALSE){
	###########################################################################
	# This function evaluates the banana function and its gradient.
	###########################################################################
	##### Inputs:
	# x:    2 x 1
	###########################################################################
	##### Outputs:
	# res$f:    			scalar
	# res$grad:				2 x 1
	# res$hess:				2 x 2
	###########################################################################
	
	x1 <- x[1];
	x2 <- x[2];
	
	res <- list();
	res$f <- 100*(x2-x1^2)^2 + (1-x1)^2;
	
	if(returnGrad){
		grad <- c(0,0);
		grad[1] <- -400*x1*(x2-x1^2) - 2*(1-x1);
		grad[2] <- 200*(x2-x1^2);
		res$grad <- grad;
	}
	
	if(returnHess){
		hess <- matrix(0,2,2);
		hess[1,1] <- -400*(x2 - 3*x1^2) + 2;
		hess[1,2] <- -400*x1;
		hess[2,1] <- hess[1,2];
		hess[2,2] <- 200;
		res$hessian <- hess;
	}
	res;
}


### NELDER-MEAD
res <- NelderMead(banana_func, x0);
res$xstar
res$fval

### NEWTON-RAPHSON
res2 <- Newton_Raphson(banana_func, x0, 1e-10, 500);
res2$xstar
res2$fval

#### BFGS
res3 <- BFGS(banana_func, x0, 1e-10, 2000);
res3$xstar
res3$fval
