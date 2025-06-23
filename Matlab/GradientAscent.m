function [xstar, fval] = GradientAscent(myfunc, x0, tolX, iterMax, stepSize, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function applies the Gradient Ascent method to maximize a function.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% myfunc:				function:
	%							takes		NumX x 1 as input
	%							returns:	objective:	scalar
	%										gradient:	NumX x 1
	% x0:					NumX x 1
	% tolX:					scalar
	% iterMax:				integer
	% stepSize:				scalar
	% varargin{1}=verbose:	boolean
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% x:					NumX x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	% Read optional arguments
	if length(varargin) >= 1; verbose  = varargin{1}; else; verbose = false; end
	
	% Do initializations
	algoName  = 'Gradient Ascent';
	x         = x0;
	iter      = 0;
	funcEvals = 0;

	while true
		iter = iter + 1;
		
		% Evaluate the function and its gradient
		[obj, grad] = myfunc(x);
		
		%%% Determine update for next step
		nextStep = stepSize * grad;
		
		%%% Compute criterionX for convergence
		criterionX = max(abs(nextStep));
		
		% Display latest results
		if verbose
			disp(sprintf('%s iteration %d: obj=%g, criterion=%g', algoName, iter, obj, criterionX));
		end

		% Output result when convergence is reached
		if criterionX <= tolX
			disp(sprintf('%s algorithm has finished afer %d iterations.', algoName, iter));
			disp(sprintf('Value of objective: %d', obj));
			disp(sprintf('Value of criterion: %d', criterionX));
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
