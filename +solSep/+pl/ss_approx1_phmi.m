function [phmi, pmd]= find_ss_approx1_phmi(sigma_0, sigma_ss, p_subset, T, pl)
%*************************************************************************
%*     Copyright c 2018 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
nsubsets = size(sigma_ss,1);
rho = sigma_ss/sigma_0;
%gamma = (pl - T)./sqrt(1+rho.^2);
gamma = (pl - T)./sqrt(sigma_ss.^2+sigma_0^2);

pmd = zeros(nsubsets,1);
for i=1:nsubsets    
    pmd(i) = normcdf(-gamma(i));
end

phmi = p_subset'*pmd;


