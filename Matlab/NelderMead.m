function [xstar, fval] = NelderMead(myfunc, x0)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function implements the Nelder-Mead method to minimize a function myfunc, taking x0 as the starting point.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% myfunc:			function (takes Nx1 as input, returns scalar)
	% x0:				N x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% xstar:			N x 1
	% fval:				scalar
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Read dimension
	N = length(x0);
	
	% Set "tuning" parameters of algorithm
	alpha = 1;
	gamma = 2;
	rho = 0.5;
	sigma = 0.5;
	
	% Set tolerance
	MaxFunEvals = 200*N;
	MaxIter     = 200*N;
	TolX   = 1e-6;
	TolFun = 1e-6;
	
	
	% Create data structure for xs and fs
	xs = zeros(N,N+1); % N x (N+1)
	fs = zeros(1,N+1); % 1 x (N+1)
	
	% Initialize values of xs
	xs(:,1) = x0;
	for ii = 1:N
		xs(:,ii+1) = x0;
%		xs(ii,ii+1) = x0(ii) + 0.01;
		if x0(ii) == 0
			xs(ii,ii+1) = 0.00025;
		else	
			xs(ii,ii+1) = 1.05*x0(ii);
		end
	end
	
	% Compute f values corresponding to initial xs
	for ii = 1:N+1
		fs(ii) = myfunc(xs(:,ii));
	end
	iter      = 0;
	funcEvals = N+1;
	fmin      = min(fs);
	
	disp(' Iteration   Func-count     min f(x)         Procedure');
	disp(sprintf('%6d %12d%17g         %s', 0, 1,fs(1), ''));
	
	message = 'initial simplex';
	
	while true
		keepInnerLoop = true;
		while keepInnerLoop
		
			%%% 1) Sort vertices by increasing values
			iter = iter + 1;
			[fs,sortidxes] = sort(fs);
			xs = xs(:,sortidxes);
			
			disp(sprintf('%6d %12d%17g         %s', iter, funcEvals,fs(1), message));

			% Check criterion and return if needed
			x_conv_crit = max(sqrt(sum((xs - xs(:,1)).^2, 1)));
			f_conv_crit = max(fs) - min(fs);

			if iter > MaxIter || funcEvals > MaxFunEvals  || (x_conv_crit < TolX && f_conv_crit < TolFun)
				xstar = xs(:,1);
				fval = fs(1);
				return;
			end
			
			
			%%% 2) Calculate centroid
			xo = mean(xs(:,1:N), 2); % N x 1
			
			%%% 3) Reflection
			xr = xo + alpha*(xo - xs(:,N+1));
			fr = myfunc(xr);
			funcEvals = funcEvals+1;
			if fr >= fs(1) && fr < fs(N)
				xs(:,N+1) = xr;
				fs(N+1) = fr;
				message = 'reflect';
				continue; % Go to 1)
			end
			
			%%% 4) Expansion
			if fr < fs(1)
				xe = xo + gamma*(xr - xo);
				fe = myfunc(xe);
				funcEvals = funcEvals+1;
				if fe < fr
					xs(:,N+1) = xe;
					fs(N+1) = fe;
					message = 'expand';
				else
					xs(:,N+1) = xr;
					fs(N+1) = fr;
					message = 'reflect';
				end
				continue; % Go to 1)
			end
			
			%%% 5) Contraction
			if fr < fs(N+1)
				xc = xo + rho*(xr - xo);
				fc = myfunc(xc);
				funcEvals = funcEvals+1;
				if fc < fr
					xs(:,N+1) = xc;
					fs(N+1) = fc;
					message = 'contract outside';
					continue; % Go to 1)
				end
			else
				xc = xo + rho*(xs(:,N+1) - xo);
				fc = myfunc(xc);
				funcEvals = funcEvals+1;
				if fc < fs(N+1)
					xs(:,N+1) = xc;
					fs(N+1) = fc;
					message = 'contract inside';
					continue; % Go to 1)
				end
			end
			keepInnerLoop = false;
		end
		
		%%% 6) Shrink
		xs(:,2:end) = xs(:,1) + sigma*(xs(:,2:end) - xs(:,1));
		for ii = 2:N+1
			fs(ii) = myfunc(xs(:,ii));
			funcEvals = funcEvals+1;
		end
		message = 'shrink';
	end	
	
	xstar = xs(:,1);
	fval = fs(1);
end
