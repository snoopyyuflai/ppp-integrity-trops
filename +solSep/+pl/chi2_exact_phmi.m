function phmi = chi2_exact_phmi(sigma_0, sigma_ss, p_subset, Kchi2, pl, dof)
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
L = pl/sigma_0;

pmd = zeros(nsubsets,1);

pmd(1) = 2*normcdf(-L);

for i=2:nsubsets    
    pmd(i) = chi2_exact_pmd(rho(i), L, Kchi2,dof);
end

phmi = p_subset'*pmd;


