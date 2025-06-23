%%%%%%%%%% PROBLEM: minimize banana function, starting from x0 as defined below %%%%%%%%%%
x0 = [-1.9; 2];

%%% NELDER-MEAD
[xstar1, fval1] = NelderMead(@banana_func, x0);
xstar1
fval1


%%% NEWTON-RAPHSON
[xstar2, fval2] = Newton_Raphson(@banana_func, x0, 1e-10, 500); % TO DO: this is for maximization, not minimization...
xstar2
fval2


%%% GRADIENT-ASCENT (with fixed step size)
[xstar2, fval2] = GradientDescent(@banana_func, x0, 1e-10, 50000, 1e-3);
xstar2
fval2


%%% BFGS
[xstar3, fval3] = BFGS(@banana_func, x0, 1e-10, 2000, true); % TO DO: this is for minimization...
xstar3
fval3


%%%%%%%%%%%%% USING CANNED CODE %%%%%%%%%%%%%
%%% Version from fminsearch (applies Nelder-Mead algorithm)
options = optimset('Display','iter','PlotFcns',@optimplotfval);
[xstar0,fval0,exitflag,output] = fminsearch(@banana_func,x0,options);
xstar0
fval0

%%% Version from fminunc (applies BFGS algorithm)
options = optimoptions(@fminunc,'SpecifyObjectiveGradient',true, ...
			'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-10, 'MaxIterations', 5e3, 'Display', 'iter');
[xstar0b, fval0b, exitflag] = fminunc(@banana_func, x0, options);
xstar0b
fval0b