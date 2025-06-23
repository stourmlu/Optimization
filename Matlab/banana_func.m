function [f, varargout] = banana_func(x)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function evaluates the banana function and its gradient.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% x:    2 x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% f:    				scalar
	% varargout{1}=grad:	2 x 1
	% varargout{2}=hess:	2 x 2
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%	disp(x');
	
	x1 = x(1);
	x2 = x(2);
	
	f = 100*(x2-x1^2)^2 + (1-x1)^2;
	
	
	if nargout >= 2
		grad = [0;0];
		grad(1) = -400*x1*(x2-x1^2) - 2*(1-x1);
		grad(2) = 200*(x2-x1^2);
		varargout{1} = grad;
	end
	
	if nargout >= 3
		hess = zeros(2,2);
		
		hess(1,1) = -400*(x2 - 3*x1^2) + 2;
		hess(1,2) = -400*x1;
		hess(2,1) = hess(1,2);
		hess(2,2) = 200;
		varargout{2} = hess;
	end

end
