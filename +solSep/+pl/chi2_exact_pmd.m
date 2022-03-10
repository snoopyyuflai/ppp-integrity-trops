function phmi = chi2_exact_pmd(rho, L, Kchi2,dof)
%*************************************************************************
%*     Copyright c 2018 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%determine search range
%value for beta = sqrt(Kchi2);

pexceed_0 = normcdf(rho*sqrt(Kchi2)-L)+normcdf(-rho*sqrt(Kchi2)-L);
pmd_chi2_0 =  ncx2cdf(Kchi2,dof,Kchi2);

% pexceed_0 = normcdf(-L)+normcdf(-L);
% pmd_chi2_0 =  ncx2cdf(Kchi2,dof,0);

%phmi_0 = pexceed_0*pmd_chi2_0;
beta_hi = (sqrt(Kchi2)+max(0,-norminv(pexceed_0)));

beta = 0:beta_hi/50:beta_hi;

pexceed = normcdf(rho*beta-L)+normcdf(-rho*beta-L);
pmd_chi2 = ncx2cdf(Kchi2,dof,beta.^2);
phmi = max(pexceed.*pmd_chi2);



% figure
% semilogy(beta, pexceed.*pmd_chi2)
% xlabel('lambda')
% ylabel('phmi')
