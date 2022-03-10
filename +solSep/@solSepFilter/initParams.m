function PARAMS = initParams(obj)




% Related to solution separation
PARAMS = struct(...
    'subsetGroupSize', 1,...
    'Psat', 1e-5,...
    'nOutSubset',1,...
    'faultResetSubsets',false,...
    'nMaxSubset',100,...
    'pl_ss_exact_optfa',0,...
    'pl_ss_approx1_optfa',1,...
    'pl_ss_approx2_optfa',0,...
    'pl_ss_exact',0,...
    'pl_ss_approx1',0,...
    'pl_ss_approx2',0,...
    'fdeMode','reinit',...  % 'reinit' (just start over when a fault is detected) or 'background'
    'pfa',1e-5,... % probability of false alert
    'phmi',1e-7,...% maximum integrity risk
    'pl_prec',0.1); % protection level precision (meters)


% List of subset setups

% Inf indicates that each element will be treated as unique
%   x specifies a particular value
%   0 means all of these are grouped for the unique sets specified elsewhere

% Example:  navsu.internal.MeasIdGnss(Inf,Inf,0,0) 
%           Measurements will be excluded per each PRN-constellation
%           (satellite) but all frequencies and measurement subtypes (code,
%           carrier, doppler) will be lumped together. 

% Example:  navsu.internal.MeasIdGnss(Inf,1,0,0)
%           This excludes all GPS satellites

PARAMS.ssBasis = [];
% GNSS measurement - all measurements from a given satellite are grouped
PARAMS.ssBasis.measId(1) = navsu.internal.MeasIdGnss(Inf,Inf,0,0);
PARAMS.ssBasis.psat(1)   = 1e-5;

% Position measurement - xyz combined 
% PARAMS.ssBasis.measId(2) = navsu.internal.MeasIdPos(Inf,0);
% PARAMS.ssBasis.psat(2)   = 1e-4;

end
