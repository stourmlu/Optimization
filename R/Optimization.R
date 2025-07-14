NelderMead <- function(myfunc, x0){
	###########################################################################
	# This function implements the Nelder-Mead method to minimize a function myfunc, taking x0 as the starting point.
	###########################################################################
	##### Inputs:
	# myfunc:			function (takes Nx1 as input, returns scalar)
	# x0:				N x 1
	###########################################################################
	##### Outputs:
	# xstar:			N x 1
	# fval:				scalar
	###########################################################################
	
	# Read dimension
	N <- length(x0);
	
	# Set "tuning" parameters of algorithm
	alpha <- 1;
	gamma <- 2;
	rho <- 0.5;
	sigma <- 0.5;
	
	# Set tolerance
	MaxFunEvals <- 200*N;
	MaxIter     <- 200*N;
	TolX   <- 1e-6;
	TolFun <- 1e-6;
	
	
	# Create data structure for xs and fs
	xs <- matrix(0,N,N+1); # N x (N+1)
	fs <- rep(0,N+1); # 1 x (N+1)
	
	# Initialize values of xs
	xs[,1] <- x0;
	for(ii in 1:N){
		xs[,ii+1] <- x0;
		if(x0[ii] == 0){
			xs[ii,ii+1] <- 0.00025;
		} else{	
			xs[ii,ii+1] <- 1.05*x0[ii];
		}
	}
	
	# Compute f values corresponding to initial xs
	for(ii in 1:(N+1)){
		fs[ii] <- myfunc(xs[,ii])$f;
	}
	iter      <- 0;
	funcEvals <- N+1;
	fmin      <- min(fs);
	
	cat(' Iteration   Func-count     min f(x)         Procedure\n');
	cat(sprintf('%6d %12d%17g         %s\n', 0, 1,fs[1], ''));
	
	message <- 'initial simplex';
	res <- list();
	
	while(TRUE){
		keepInnerLoop <- TRUE;
		while(keepInnerLoop){
		
			### 1) Sort vertices by increasing values
			iter <- iter + 1;
			mysort <- sort(fs, index.return=TRUE);
			fs <- mysort$x;
			xs <- xs[,mysort$ix];
			
			cat(sprintf('%6d %12d%17g         %s\n', iter, funcEvals,fs[1], message));

			# Check criterion and return if needed
			x_conv_crit <- max(sqrt(apply((xs - xs[,1])^2, 2, sum)));
			f_conv_crit <- max(fs) - min(fs);

			if(iter > MaxIter || funcEvals > MaxFunEvals  || (x_conv_crit < TolX && f_conv_crit < TolFun)){
				res$xstar <- xs[,1];
				res$fval  <- fs[1];
				return(res);
			}
			
			### 2) Calculate centroid
			xo <- apply(xs[,1:N], 1, mean); # N x 1
			
			### 3) Reflection
			xr <- xo + alpha*(xo - xs[,N+1]);
			fr <- myfunc(xr)$f;
			funcEvals <- funcEvals+1;
			if(fr >= fs[1] && fr < fs[N]){
				xs[,N+1] <- xr;
				fs[N+1]  <- fr;
				message <- 'reflect';
				next; # Go to 1)
			}
			
			### 4) Expansion
			if(fr < fs[1]){
				xe <- xo + gamma*(xr - xo);
				fe <- myfunc(xe)$f;
				funcEvals <- funcEvals+1;
				if(fe < fr){
					xs[,N+1] <- xe;
					fs[N+1]  <- fe;
					message <- 'expand';
				} else{
					xs[,N+1] <- xr;
					fs[N+1]  <- fr;
					message <- 'reflect';
				}
				next; # Go to 1)
			}
			
			### 5) Contraction
			if(fr < fs[N+1]){
				xc <- xo + rho*(xr - xo);
				fc <- myfunc(xc)$f;
				funcEvals <- funcEvals+1;
				if(fc < fr){
					xs[,N+1] <- xc;
					fs[N+1]  <- fc;
					message <- 'contract outside';
					next; # Go to 1)
				}
			} else {
				xc <- xo + rho*(xs[,N+1] - xo);
				fc <- myfunc(xc)$f;
				funcEvals <- funcEvals+1;
				if(fc < fs[N+1]){
					xs[,N+1] <- xc;
					fs[N+1]  <- fc;
					message <- 'contract inside';
					next; # Go to 1)
				}
			}
			keepInnerLoop <- FALSE;
		}
		
		### 6) Shrink
		xs[,2:(N+1)] <- xs[,1] + sigma*(xs[,2:(N+1)] - xs[,1]);
		for(ii in 2:N+1){
			fs[ii] <- myfunc(xs[,ii])$f;
			funcEvals <- funcEvals+1;
		}
		message <- 'shrink';
	}	
	res$xstar <- xs[,1];
	res$fval  <- fs[1];
}


Newton_Raphson <- function(myfunc, x0, tolX, iterMax, stepSize=1, verbose=FALSE){
	###########################################################################
	# This function applies the Newton-Raphson method to optimize a function.
	###########################################################################
	##### Inputs:
	# myfunc:				function:
	#							takes		NumX x 1 as input
	#							returns:	objective:	scalar
	#										gradient:	NumX x 1
	#										hessian:	NumX x NumX
	# x0:					NumX x 1
	# tolX:					scalar
	# iterMax:				integer
	# stepSize:				scalar (between 0 and 1)
	# verbose:				boolean
	###########################################################################
	##### Outputs:
	# x:					NumX x 1
	###########################################################################
	
	# Do initializations
	algoName  <- 'Newton-Raphson';
	x         <- x0;
	iter      <- 0;
	funcEvals <- 0;

	# Iterate
	while(TRUE){
		iter <- iter + 1;
		
		# Evaluate the function, its gradient and its Hessian
		tmp <- myfunc(x, TRUE, TRUE);
		obj  <- tmp$f;
		grad <- tmp$grad;
		hess <- tmp$hess;
		
		### Determine update for next step
		nextStep <- -stepSize * solve(hess,grad);
		
		### Compute criterionX for convergence
		criterionX <- max(abs(nextStep));
		
		# Display latest results
		if(verbose){
			cat(sprintf('%s iteration %d: obj=%g, criterion=%g\n', algoName, iter, obj, criterionX));
		}
		
		# Output result when convergence is reached or iterMax is reached
		if(criterionX <= tolX || iter >= iterMax){
			if(criterionX <= tolX){
				cat(sprintf('%s algorithm has finished afer %d iterations.\n', algoName, iter));
			} else{
				cat(sprintf('%s algorithm has not converged afer %d iterations.\n', algoName, iter));
			}
			cat(sprintf('Value of objective: %g\n', obj));
			cat(sprintf('Value of criterion: %g\n', criterionX));
			res <- list();
			res$xstar <- x;
			res$fval  <- obj;
			res$grad  <- grad;
			res$hess  <- hess;
			return(res);
		}
		
		# Prepare next round
		x <- x + nextStep;
	}
}


BFGS <- function(myfunc, x0, tolX, iterMax, verbose=FALSE){
	###########################################################################
	# This function minimizes a function via the BFGS method, with a cubic interpolation to perform the line search.
	###########################################################################
	##### Inputs:
	# myfunc:				function:
	#							takes		NumX x 1 as input
	#							returns:	objective:	scalar
	#										gradient:	NumX x 1
	# x0:					NumX x 1
	# tolX:					scalar
	# iterMax:				integer
	# verbose:	boolean
	###########################################################################
	##### Outputs:
	# xstar
	# fval
	# grad
	###########################################################################
	
	# Read dimensions
	NumX <- length(x0);
	identity <- diag(NumX);

	# Do initializations
	iter <- 0;
	criterionX <- Inf;
	
	# Initialize x, grad and B_inv
	x_old <- x0;
	tmp <- myfunc(x_old, TRUE);
	obj_old  <- tmp$f;
	grad_old <- tmp$grad;
	
	
	funcCount <- 1;
	B_inv     <- diag(NumX);
	alpha     <- 1/max(abs(grad_old))
	
	# Display stuff
	if(verbose){
		disp('                                                        First-order');
		disp(' Iteration  Func-count       f(x)        Step-size       optimality');
			firstOrderOptimality <- max(abs(grad_old));
			cat(sprintf('%6d%12d%17g%30.3g\n', 0, funcCount,obj_old,firstOrderOptimality));
	}
	
	while(iter < iterMax && criterionX > tolX){
		iter <- iter + 1;
		
		# 1) Obtain direction
		direction <- -B_inv %*% grad_old;
		
		foundGoodCandidate <- FALSE;
		while(!foundGoodCandidate){
			delta_x <- as.numeric(alpha) * direction;

			if(max(abs(delta_x)) < tolX){
				x_new    <- x_old;
				obj_new  <- obj_old;
				grad_new <- grad_old;
				foundGoodCandidate <- TRUE;
				break;
			}
			
			x_new <- x_old + delta_x;
			tmp   <- myfunc(x_new, TRUE);
			obj_new  <- tmp$f;
			grad_new <- tmp$grad;
			
			delta_grad          <- grad_new - grad_old;
			funcCount           <- funcCount + 1;
			
			obj_diff            <- obj_new - obj_old;
			grad_new_TIMES_dir  <- t(grad_new) %*% direction;
			
			myalpha <- alpha; # Store alpha for Display

			### Do cubic interpolation to find alpha_c
			gprime0      <- t(grad_old) %*% direction;
			gprime_alpha <- grad_new_TIMES_dir;
			g0      <- obj_old;
			galpha  <- obj_new;
			beta1   <- gprime0 + gprime_alpha + 3*(g0 - galpha)/(alpha);
			beta2   <- sqrt(max(0,beta1^2 - gprime0*gprime_alpha)); ### Including max(0,.) avoids taking the sqrt of negative values
			alpha_c <- alpha * (1 - (gprime_alpha + beta2 - beta1)/(gprime_alpha - gprime0 + 2*beta2));
			
			if(alpha_c < 0){
				alpha_c <- 2*alpha;
			}
			
			if(obj_diff >= 0){
				if(grad_new_TIMES_dir > 0){
					# Case 1
					if(alpha < 0.1){
						alpha <- 0.5*alpha_c;
					} else{
						alpha <- alpha_c;
					}
				} else{
					# Case 4
					alpha <- min(alpha_c, 0.5*alpha);
				}
			} else{
				denom <- t(delta_x) %*% delta_grad;
				if(denom >= 0){
					# Perform update
					denom_inv <- as.numeric(1/denom);
					tmp1      <- denom_inv * (delta_x %*% t(delta_grad)); # NumX x NumX
					B_inv <- (identity - tmp1) %*% B_inv %*% (identity - t(tmp1)) + denom_inv * delta_x %*% t(delta_x);
					if(grad_new_TIMES_dir >= 0){
						# Case 2a
						alpha <- min(1, alpha_c);
					} else{
						# Case 3a
						p <- min(1 + denom - grad_new_TIMES_dir + min(0,alpha));
						alpha <- min(min(2,p),1.2*alpha_c);
					}
					foundGoodCandidate <- TRUE;
					
				} else{
					if(grad_new_TIMES_dir >= 0){
						# Case 2b
						alpha <- 0.9*alpha_c;
					} else{
						# Case 3b
						alpha <- min(2, max(1.5, alpha));
						alpha <- min(alpha, alpha_c);
					}
				}
			}
		}
		
		# Update x, grad, criterion
		x_old      <- x_new;
		obj_old    <- obj_new;
		grad_old   <- grad_new;
		criterionX <- max(abs(delta_x));
		
		# Display stuff
		if(verbose){
			firstOrderOptimality <- max(abs(grad_new));
			disp(sprintf('%6d%12d%17g%15g%15.3g', iter, funcCount,obj_new,myalpha,firstOrderOptimality));
		}	
	}
	cat(sprintf('BFGS algorithm has finished afer %d iterations.\n', iter));
	cat(sprintf('Value of objective: %g\n', obj_new));
	cat(sprintf('Value of criterion: %g\n', criterionX));
	
	# Output everything
	res <- list();
	res$xstar <- x_new;
	res$fval  <- obj_new;
	res$grad  <- grad_new;
	return(res);
}
