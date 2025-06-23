function [xstar, fval] = Newton_Raphson(myfunc, x0, tolX, iterMax, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function applies the Newton-Raphson method to optimize a function.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% myfunc:				function:
	%							takes		NumX x 1 as input
	%							returns:	objective:	scalar
	%										gradient:	NumX x 1
	%										hessian:	NumX x NumX
	% x0:					NumX x 1
	% tolX:					scalar
	% iterMax:				integer
	% varargin{1}=stepSize:	scalar (between 0 and 1)
	% varargin{2}=verbose:	boolean
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% x:					NumX x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Read optional arguments
	if length(varargin) >= 1; stepSize = varargin{1}; else; stepSize = 1;    end
	if length(varargin) >= 2; verbose  = varargin{2}; else; verbose = false; end
	
	% Do initializations
	algoName  = 'Newton-Raphson';
	x         = x0;
	iter      = 0;
	funcEvals = 0;

	% Iterate
	while true
		iter = iter + 1;
		
		% Evaluate the function, its gradient and its Hessian
		[obj, grad, hess] = myfunc(x);
		
		%%% Determine update for next step
		nextStep = -stepSize * hess\grad;
		
		%%% Compute criterionX for convergence
		criterionX = max(abs(nextStep));
		
		% Display latest results
		if verbose
			disp(sprintf('%s iteration %d: obj=%g, criterion=%g', algoName, iter, obj, criterionX));
		end
		
		% Output result when convergence is reached
		if criterionX <= tolX
			disp(sprintf('%s algorithm has finished afer %d iterations.', algoName, iter));
			disp(sprintf('Value of objective: %g', obj));
			disp(sprintf('Value of criterion: %g', criterionX));
			xstar = x;
			fval  = obj;
			return;
		end
		
		% Output result when iterMax is reached
		if iter >= iterMax
			disp(sprintf('%s algorithm has not converged afer %d iterations.', algoName, iter));
			xstar = x;
			fval = obj;
			return;
		end
		
		% Prepare next round
		x = x + nextStep;
	end
end
