function [x, fval, exitflag] = numjac_Newton_method(fun, x0, delta_numjac, ...
    alpha, ...
    max_fun_evals, ...
    tolx, ...
    tolf)
%
% A humble N-dimensional Newton method for when we encounter the 'Maximum
% number of users for Optimization_Toolbox reached' error... Solves
% fun(x) = 0, with x being an N-d vector, with damping factor alpha (multiplies
% the step) and seed x0. Calculates the numerical jacobian with a perturbation
% size equal to the input delta_numjac. A typical value for it is sqrt(eps).
% Returns [x, fval, exitflag].
%
% exitflag can take the following values
%     * 1 (numeric_success) -> Function converged to a solution x with abs(fun(x)) <= tolf
%     * -2                  -> Function converged to a solution x with abs(fun(x)) > tolf
%     * -3                  -> NaN or Inf function value was encountered
%     * -100                -> limit of max. numer of function evaluations was reached
%

exitflag = 1;
x        = x0;
fval     = fun(x);

if ~all(isfinite(x)) || ~all(isfinite(fval))
    exitflag = -3;
    x        = nan;
    fval     = nan;
    
elseif sumabs(fval) > tolf
    n = numel(x0);
    
    for fun_evals = 2:(2 * n):max_fun_evals
        
        numjac = zeros(n);
        for i = 1:n
            x(i)         = x(i) + delta_numjac;
            f_p          = fun(x);
            x(i)         = x(i) - 2 * delta_numjac;
            f_m          = fun(x);
            numjac(:, i) = 0.5 * (f_p - f_m) / delta_numjac;
            x(i)         = x(i) + delta_numjac;

        end
        
        dx   = numjac \ fval;
        x    = x - alpha * dx;
        fval = fun(x);
        
        if sumabs(fval) < tolf
            return

        elseif sumabs(0.5 * dx) <= 2.0 * tolx * max(1.0, sumabs(x))
            exitflag = -2;
            return
            
        elseif ~all(isfinite(x)) || ~all(isfinite(fval))
            exitflag = -3;
            x        = nan;
            fval     = nan;
            return
            
        end
        
    end
    
    exitflag = -100;
    
end

end