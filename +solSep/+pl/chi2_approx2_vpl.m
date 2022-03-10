function pl = chi2_approx2_vpl(sigma_0, sigma_ss, p_subset, Kchi2, phmi)
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

K_hi = - norminv(phmi./(p_subset*nsubsets));
pl_hi = sqrt(Kchi2)*max(sigma_ss) +max(K_hi.*sqrt(sigma_0^2+sigma_ss.^2));

K0 = - norminv(.5*phmi/nsubsets);
pl_lo = K0*sigma_0;

pl = max(pl_hi, pl_lo);





