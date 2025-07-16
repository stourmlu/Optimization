include("Optimization.jl")

########## PROBLEM: minimize banana function, starting from x0 as defined below ##########
x0 = [-1.9; 2]

function banana_func(x; grad::Bool=false, hessian::Bool=false)
    ###########################################################################
    # This function evaluates the banana function and its gradient.
    ###########################################################################
    ##### Inputs:
    # x:    2 x 1
    ###########################################################################
    ##### Outputs:
    # f:    				scalar
    ###########################################################################

    x1 = x[1]
    x2 = x[2]
    
    f = 100*(x2-x1^2)^2 + (1-x1)^2
	
	if grad
		g = zeros(2)
		g[1] = -400*x1*(x2-x1^2) - 2*(1-x1)
		g[2] = 200*(x2-x1^2)
	else
		g = nothing
	end
	
	if hessian
		hess = zeros(2,2)
		hess[1,1] = -400*(x2 - 3*x1^2) + 2
		hess[1,2] = -400*x1
		hess[2,1] = hess[1,2]
		hess[2,2] = 200
	else
		hess = nothing
	end
	
	return (f, g, hess)
end

### NELDER-MEAD
(xstar1, fval1) = NelderMead(x0, banana_func)

### NEWTON-RAPHSON
(xstar2, fval2) = Newton_Raphson(x0, banana_func)

### BFGS
(xstar3, fval3) = BFGS(x0, banana_func)
