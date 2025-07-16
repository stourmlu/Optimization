using Printf
using Statistics
using LinearAlgebra

function NelderMead(x0, myfunc)
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
	N = length(x0)
	
	# Set "tuning" parameters of algorithm
	alpha = 1
	gamma = 2
	rho = 0.5
	sigma = 0.5
	
	# Set tolerance
	MaxFunEvals = 200*N
	MaxIter     = 200*N
	TolX   = 1e-6
	TolFun = 1e-6
	
	
	# Create data structure for xs and fs
	xs = zeros(N,N+1) # N x (N+1)
	fs = zeros(N+1) # 1 x (N+1)
	
	# Initialize values of xs
	xs[:,1] = x0
	for ii = 1:N
		xs[:,ii+1] = x0
		if x0[ii] == 0
			xs[ii,ii+1] = 0.00025
		else
			xs[ii,ii+1] = 1.05*x0[ii]
		end
	end
	
	# Compute f values corresponding to initial xs
	for ii = 1:N+1
		fs[ii] = myfunc(xs[:,ii])
	end
	iter      = 0
	funcEvals = N+1
	fmin      = minimum(fs)
	
	println(" Iteration   Func-count     min f(x)         Procedure")
	println(@sprintf("%6d %12d%17g         %s", 0, 1,fs[1], ""))
	
	message = "initial simplex"
	
	sortidxes = zeros(Int, N+1)
	
	while true
		keepInnerLoop = true
		while keepInnerLoop
		
			### 1) Sort vertices by increasing values
			iter = iter + 1
			sortidxes = sortperm!(sortidxes, fs)
			fs = fs[sortidxes]
			xs = xs[:,sortidxes]
			
			println(@sprintf("%6d %12d%17g         %s", iter, funcEvals,fs[1], message))

			# Check criterion and return if needed
			x_conv_crit = maximum(sqrt.(sum((xs .- xs[:,1]).^2, dims=2)))
			f_conv_crit = maximum(fs) - minimum(fs)

			if iter > MaxIter || funcEvals > MaxFunEvals  || (x_conv_crit < TolX && f_conv_crit < TolFun)
				xstar = xs[:,1]
				fval = fs[1]
				return xstar, fval
			end
			
			
			### 2) Calculate centroid
			xo = Statistics.mean(xs[:,1:N], dims=2) # N x 1
			
			### 3) Reflection
			xr = xo + alpha*(xo - xs[:,N+1])
			fr = myfunc(xr)
			funcEvals = funcEvals+1
			if fr >= fs[1] && fr < fs[N]
				xs[:,N+1] = xr
				fs[N+1] = fr
				message = "reflect"
				continue # Go to 1)
			end
			
			### 4) Expansion
			if fr < fs[1]
				xe = xo + gamma*(xr - xo)
				fe = myfunc(xe)
				funcEvals = funcEvals+1
				if fe < fr
					xs[:,N+1] = xe
					fs[N+1] = fe
					message = "expand"
				else
					xs[:,N+1] = xr
					fs[N+1] = fr
					message = "reflect"
				end
				continue # Go to 1)
			end
			
			### 5) Contraction
			if fr < fs[N+1]
				xc = xo + rho*(xr - xo)
				fc = myfunc(xc)
				funcEvals = funcEvals+1
				if fc < fr
					xs[:,N+1] = xc
					fs[N+1] = fc
					message = "contract outside"
					continue # Go to 1)
				end
			else
				xc = xo + rho*(xs[:,N+1] .- xo)
				fc = myfunc(xc)
				funcEvals = funcEvals+1
				if fc < fs[N+1]
					xs[:,N+1] = xc
					fs[N+1] = fc
					message = "contract inside"
					continue # Go to 1)
				end
			end
			keepInnerLoop = false
		end
		
		### 6) Shrink
		xs[:,2:end] = xs[:,1] .+ sigma*(xs[:,2:end] .- xs[:,1])
		for ii = 2:N+1
			fs[ii] = myfunc(xs[:,ii])
			funcEvals = funcEvals+1
		end
		message = "shrink"
	end
	
	xstar = xs[:,1]
	fval = fs[1]
	return xstar, fval
end



function Newton_Raphson(x0, myfunc; tolX=1e-10, iterMax=500, stepSize=1, verbose=false)
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
	algoName  = "Newton-Raphson"
	x         = x0
	iter      = 0
	funcEvals = 0

	# Iterate
	while true
		iter = iter + 1
		
		# Evaluate the function, its gradient and its Hessian
		(obj, grad, hess) = myfunc(x, grad=true, hessian=true)
		
		### Determine update for next step
		nextStep = -stepSize * hess\grad
		
		### Compute criterionX for convergence
		criterionX = maximum(abs.(nextStep))
		
		# Display latest results
		if verbose
			println(@sprintf("%s iteration %d: obj=%g, criterion=%g", algoName, iter, obj, criterionX))
		end
		
		# Output result when convergence is reached or iterMax is reached
		if criterionX <= tolX || iter >= iterMax
			if criterionX <= tolX
				println(@sprintf("%s algorithm has finished afer %d iterations.", algoName, iter))
			else
				println(@sprintf("%s algorithm has not converged afer %d iterations.", algoName, iter))
			end
			println(@sprintf("Value of objective: %g", obj))
			println(@sprintf("Value of criterion: %g", criterionX))
			xstar = x
			fval  = obj
			return xstar, fval, grad, hess
		end
		
		# Prepare next round
		x = x + nextStep
	end
end


function BFGS(x0, myfunc; tolX=1e-10, iterMax=500, stepSize=1, verbose=false)
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
	# verbose:				boolean
	###########################################################################
	##### Outputs:
	# xstar
	# fval
	# grad
	###########################################################################
	
	# Read dimensions
	NumX = length(x0)

	# Do initializations
	iter = 0
	criterionX = Inf
	
	# Initialize x, grad and B_inv
	x_old = x0
	(obj_old,grad_old) = myfunc(x_old, grad=true)
	funcCount = 1
	B_inv = I
	alpha = 1. /maximum(abs.(grad_old))
	
	# Display stuff
	if verbose
	println("                                                        First-order")
	println(" Iteration  Func-count       f(x)        Step-size       optimality")
		firstOrderOptimality = maximum(abs.(grad_old))
		println(@sprintf("%6d%12d%17g%30.3g", 0, funcCount,obj_old,firstOrderOptimality))
	end

	##### Allocate memory for what goes in the loop #####
	x_new = x0
	obj_new = obj_old
	grad_new = grad_old
	delta_x = grad_old
	#####################################################
	while iter < iterMax && criterionX > tolX
		iter = iter + 1
		
		# 1) Obtain direction
		direction = -B_inv * grad_old
		
		foundGoodCandidate = false
		while !foundGoodCandidate
			delta_x = alpha * direction
			
			if maximum(abs.(delta_x)) < tolX
				x_new    = x_old
				obj_new  = obj_old
				grad_new = grad_old
				foundGoodCandidate = true
				break
			end
			
			x_new               = x_old + delta_x
			(obj_new, grad_new) = myfunc(x_new, grad=true)
			delta_grad          = grad_new - grad_old
			funcCount           = funcCount + 1
			
			obj_diff            = obj_new - obj_old
			grad_new_TIMES_dir  = grad_new'*direction

			### Do cubic interpolation to find alpha_c
			gprime0      = grad_old'*direction
			gprime_alpha = grad_new_TIMES_dir
			g0      = obj_old
			galpha  = obj_new
			beta1   = gprime0 + gprime_alpha + 3*(g0 - galpha)/(alpha)
			beta2   = sqrt(max(0,beta1^2 - gprime0*gprime_alpha)) ### Including max(0,.) avoids taking the sqrt of negative values
			alpha_c = alpha * (1 - (gprime_alpha + beta2 - beta1)/(gprime_alpha - gprime0 + 2*beta2))

			if alpha_c < 0
				alpha_c = 2*alpha
			end
			
			if obj_diff >= 0
				if grad_new_TIMES_dir > 0
					# Case 1
					if alpha < 0.1
						alpha = 0.5*alpha_c
					else
						alpha = alpha_c
					end
				else
					# Case 4
					alpha = min(alpha_c, 0.5*alpha)
				end
			else
				denom = delta_x'*delta_grad
				if denom >= 0
					# Perform update
					denom_inv = 1/denom
					tmp1      = denom_inv * (delta_x*delta_grad') # NumX x NumX
					B_inv = (I - tmp1) * B_inv * (I - tmp1') + denom_inv * delta_x*delta_x'
					if grad_new_TIMES_dir >= 0
						# Case 2a
						alpha = min(1, alpha_c)
					else
						# Case 3a
						p = min(1 + denom - grad_new_TIMES_dir + min(0,alpha))
						alpha = min(min(2,p),1.2*alpha_c)
					end
					foundGoodCandidate = true
					
				else
					if grad_new_TIMES_dir >= 0
						# Case 2b
						alpha = 0.9*alpha_c
					else
						# Case 3b
						alpha = min(2, max(1.5, alpha))
						alpha = min(alpha, alpha_c)
					end
				end
			end
		end
		
		# Update x, grad, criterion
		x_old      = x_new
		obj_old    = obj_new
		grad_old   = grad_new
		criterionX = maximum(abs.(delta_x))
		
		# Display stuff
		if verbose
			firstOrderOptimality = maximum(abs.(grad_new))
			println(@sprintf("%6d%12d%17g%15g%15.3g", iter, funcCount,obj_new,alpha,firstOrderOptimality))
		end
		
	end
	println(@sprintf("BFGS algorithm has finished afer %d iterations.", iter))
	println(@sprintf("Value of objective: %g", obj_new))
	println(@sprintf("Value of criterion: %g", criterionX))
	
	# Output everything
	xstar = x_new
	fval = obj_new
	grad = grad_new
	
	return xstar, fval, grad
end
