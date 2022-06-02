function CAS_mps = cas(Mach, pressure_Pa, density_SL_kgpm3, pressure_SL_Pa)
%
% Returns the calibrated airspeed (CAS) given the true Mach number (formed
% with the true airspeed TAS) and the atmospheric pressure. User can
% optionally specify the sea level air density and pressure. All in SI
% units.
%
% Algorithm:
%
%  * The CAS is the speed that the Pitot tube measures. A Pitot tube
%    obtains a speed from the dynamic_pressure == impact_pressure =
%    (stagnation_pressure - static_pressure).
%
%  * Isentropic flow (with Rankine-Hugoniot for normal shockwaves in
%    supersonic regime) formulas yield
%
%         impact_pressure = F(Mach, pressure),
%
%    and the CAS, by definition, can be obtained by solving
%
%         impact_pressure = F(CAS / soundspeed_SL, pressure_SL).
%

% validate inputs
no_density_SL  = nargin < 3 || isempty(density_SL_kgpm3);
no_pressure_SL = nargin < 4 || isempty(pressure_SL_Pa);
if no_density_SL
    if no_pressure_SL
        [density_SL_kgpm3, ~, ~, ~, pressure_SL_Pa] = isa_atmos();
    else
        density_SL_kgpm3 = isa_atmos();
    end
    
end

if ~no_density_SL && no_pressure_SL
    [~, ~, ~, ~, pressure_SL_Pa] = isa_atmos();
end


% constants
gamma = 1.4;
gm1   = gamma - 1.;
gp1   = gamma + 1.;

soundspeed2_SL = gamma * pressure_SL_Pa / density_SL_kgpm3;
Mach2          = Mach * Mach;


% sonic impact pressure evaluated at SL and Mach = 1.
qlimit_SL = pressure_SL_Pa * ((1. + 0.5 * gm1)^(gamma / gm1) - 1.);


% evaluate impact pressure
if Mach <= 1.0
    impact_pressure = pressure_Pa * ((1. + 0.5 * gm1 * Mach2)^(gamma / gm1) - 1.);
    
else
    p_quot_Rankine    = 1. + 2. * gamma / gp1 * (Mach2 - 1.);                                % p_postshockwave / p_flight
    p_quot_isentropic = (gp1 * gp1 * Mach2 / (4. * gamma * Mach2 - 2. * gm1))^(gamma / gm1); % p_stagnation / p_postshockwave
    
    impact_pressure = pressure_Pa * (p_quot_isentropic * p_quot_Rankine - 1.);
    
end


% evaluate CAS
if impact_pressure <= qlimit_SL
    CAS_mps = sqrt(2. * soundspeed2_SL / gm1 * ...
        ((impact_pressure / pressure_SL_Pa + 1.)^(gm1 / gamma) - 1.));
    
else
    % the following equation must be solved:
    % k * d = y / (1 - b / y)^n
    % y = calibrated_Mach^2 = (cas / soundspeed_SL)^2
    % d = impact_pressure / pressure_SL + 1
    % k, b, n = f(gamma)
    
    k = 1. / (0.5 * gp1)^(gamma / gm1) / (0.5 * gp1 / gamma)^(1 / gm1); % 1. / (1.2^3.5 * (6 / 7)^2.5)
    b = 0.5 * gm1 / gamma;
    n = 1. / gm1;
    
    d    = impact_pressure / pressure_SL_Pa + 1.;
    kd   = k * d;
    kdnb = kd * n * b;
    
    max_evals            = 30;
    tol                  = 1e-8;
    calibrated_Mach_seed = 1.0;
    seed                 = calibrated_Mach_seed * calibrated_Mach_seed;
    
    fun     = @(y)(kd * (1 - b / y)^n - y);
    dfun_dy = @(y)(kdnb * (1. - b / y)^(n - 1.) / y / y - 1.);
    
    [y, ~, exitflag] = Newton_method(fun, dfun_dy, seed, 1., ...
        max_evals, ...
        tol, ...
        tol);
    
    if exitflag ~= 1
        warning('cas: the root finding algorithm failed, result will be approximate.');
        
        % regression (I found it in the EFA code...)
        y = 0.41726 + 0.7767 * (d - 1.) - 0.0989 / (d - 1.);
        
    end
    
    CAS_mps = sqrt(soundspeed2_SL * y);
    
end

end