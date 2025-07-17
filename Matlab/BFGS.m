function [xstar,fval,grad] = BFGS(myfunc, x0, tolX, iterMax, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function minimizes a function via the BFGS method, with a cubic interpolation to perform the line search.
	%
	% Written based on the explanations available here:
	%	http://www.ece.northwestern.edu/local-apps/matlabhelp/toolbox/optim/fminunc.html
	%	http://www.ece.northwestern.edu/local-apps/matlabhelp/toolbox/optim/tutori6b.html#166
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% myfunc:				function:
	%							takes		NumX x 1 as input
	%							returns:	objective:	scalar
	%										gradient:	NumX x 1
	% x0:					NumX x 1
	% tolX:					scalar
	% iterMax:				integer
	% varargin{1}=verbose:	boolean
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% xstar
	% fval
	% grad
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	% Read optional arguments
	if length(varargin) >= 1; verbose  = varargin{1}; else; verbose = false; end
	
	% Read dimensions
	NumX = length(x0);
	identity = eye(NumX);

	% Do initializations
	iter = 0;
	criterionX = Inf;
	
	% Initialize x, grad and B_inv
	x_old = x0;
	[obj_old,grad_old] = myfunc(x_old);
	funcCount = 1;
	B_inv = eye(NumX);
	alpha = 1./max(abs(grad_old));
	
	% Display stuff
	if verbose
	disp('                                                        First-order');
	disp(' Iteration  Func-count       f(x)        Step-size       optimality');
		firstOrderOptimality = max(abs(grad_old));
		disp(sprintf('%6d%12d%17g%30.3g', 0, funcCount,obj_old,firstOrderOptimality));
	end
	
	while iter < iterMax && criterionX > tolX
		iter = iter + 1;
		
		% 1) Obtain direction
		direction = -B_inv * grad_old;
		
		foundGoodCandidate = false;
		while ~foundGoodCandidate
			delta_x = alpha * direction;

			if max(abs(delta_x)) < tolX
				x_new    = x_old;
				obj_new  = obj_old;
				grad_new = grad_old;
				foundGoodCandidate = true;
				break;
			end
			
			x_new               = x_old + delta_x;
			[obj_new, grad_new] = myfunc(x_new);
			delta_grad          = grad_new - grad_old;
			funcCount           = funcCount + 1;
			
			obj_diff            = obj_new - obj_old;
			grad_new_TIMES_dir  = grad_new'*direction;
			
			myalpha = alpha; % Store alpha for display

			%%% Do cubic interpolation to find alpha_c
			gprime0      = grad_old'*direction;
			gprime_alpha = grad_new_TIMES_dir;
			g0      = obj_old;
			galpha  = obj_new;
			beta1   = gprime0 + gprime_alpha + 3*(g0 - galpha)/(alpha);
			beta2   = sqrt(max(0,beta1^2 - gprime0*gprime_alpha)); %%% Including max(0,.) avoids taking the sqrt of negative values
			alpha_c = alpha * (1 - (gprime_alpha + beta2 - beta1)/(gprime_alpha - gprime0 + 2*beta2));

			if alpha_c < 0
				alpha_c = 2*alpha;
			end
			
			if obj_diff >= 0
				if grad_new_TIMES_dir > 0
					% Case 1
					if alpha < 0.1
						alpha = 0.5*alpha_c;
					else
						alpha = alpha_c;
					end
				else
					% Case 4
					alpha = min(alpha_c, 0.5*alpha);
				end
			else
				denom = delta_x'*delta_grad;
				if denom >= 0
					% Perform update
					denom_inv = 1/denom;
					tmp1      = denom_inv * (delta_x*delta_grad'); % NumX x NumX
					B_inv = (identity - tmp1) * B_inv * (identity - tmp1') + denom_inv * delta_x*delta_x';
					if grad_new_TIMES_dir >= 0
						% Case 2a
						alpha = min(1, alpha_c);
					else
						% Case 3a
						p = min(1 + denom - grad_new_TIMES_dir + min(0,alpha));
						alpha = min(min(2,p),1.2*alpha_c);
					end
					foundGoodCandidate = true;
					
				else
					if grad_new_TIMES_dir >= 0
						% Case 2b
						alpha = 0.9*alpha_c;
					else
						% Case 3b
%						direction = -grad_old; % "Change to steepest descent momentarily" --> ?
						alpha = min(2, max(1.5, alpha));
						alpha = min(alpha, alpha_c);
					end
				end
			end
		end
		
		% Update x, grad, criterion
		x_old      = x_new;
		obj_old    = obj_new;
		grad_old   = grad_new;
		criterionX = max(abs(delta_x));
		
		% Display stuff
		if verbose
			firstOrderOptimality = max(abs(grad_new));
			disp(sprintf('%6d%12d%17g%15g%15.3g', iter, funcCount,obj_new,myalpha,firstOrderOptimality));
		end
		
	end
	disp(sprintf('BFGS algorithm has finished afer %d iterations.', iter));
	disp(sprintf('Value of objective: %g', obj_new));
	disp(sprintf('Value of criterion: %g', criterionX));
	
	% Output everything
	xstar = x_new;
	fval = obj_new;
	grad = grad_new;
end
