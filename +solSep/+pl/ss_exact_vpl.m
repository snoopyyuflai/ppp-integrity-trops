function pl_hi = ss_exact_vpl(sigma_0, sigma_ss, p_subset, T, phmi,tol)
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

%initial values:
K_hi = - norminv(.5*phmi./(p_subset*nsubsets));
pl_hi = max(T+K_hi.*sqrt(sigma_ss.^2 + sigma_0^2));

K0 = - norminv(.5*phmi);
pl_lo = K0*sigma_0;

count= 1;
%evaluate pl_hi
while (abs(pl_hi-pl_lo)>tol)&&(count<100)
   count = count+1;
   pl_mid = .5*(pl_lo+pl_hi);
   phmi_mid = solSep.pl.ss_exact_phmi(sigma_0, sigma_ss, p_subset, T, pl_mid);
   if phmi_mid>phmi
       pl_lo = pl_mid;
   else
       pl_hi = pl_mid;
   end
end





