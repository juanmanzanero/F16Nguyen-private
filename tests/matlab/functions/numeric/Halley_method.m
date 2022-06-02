function [x, fval, exitflag, num_fun_evals] = Halley_method(fun, dfun_dx, d2fun_dx2, ...
    x0, ...
    alpha, ...
    max_fun_evals, ...
    tolx, ...
    tolf)
%
% Halley's method to solve fun(x) = 0, with damping factor alpha
% (multiplies the step) and seed x0. Returns [x, fval, exitflag, fun_evals].
%
% exitflag can take the following values
%     * 1 (numeric_success) -> Function converged to a solution x with abs(fun(x)) <= tolf
%     * -2                  -> Function converged to a solution x with abs(fun(x)) > tolf
%     * -3                  -> NaN or Inf function value was encountered
%     * -100                -> limit of max. numer of function evaluations was reached
%

exitflag      = 1;
x             = x0;
fval          = fun(x);
num_fun_evals = 1;
if ~isfinite(x) || ~isfinite(fval)
    exitflag = -3;
    x        = nan;
    fval     = nan;
    
elseif abs(fval) > tolf
    for num_fun_evals = 2:max_fun_evals
        dfval = dfun_dx(x);
        dx    = 2 * fval * dfval / (2 * dfval * dfval - fval * d2fun_dx2(x));
        x     = x - alpha * dx;
        fval  = fun(x);
        
        if (abs(fval) < tolf)
            return
            
        elseif abs(0.5 * dx) <= 2.0 * tolx * max(1.0, abs(x))
            exitflag = -2;
            return
            
        elseif ~isfinite(x) || ~isfinite(fval)
            exitflag = -3;
            x        = nan;
            fval     = nan;
            return
            
        end
        
    end
    
    exitflag = -100;
    
end

end