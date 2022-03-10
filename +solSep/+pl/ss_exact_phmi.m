function [phmi, pmd]= ss_exact_phmi(sigma_0, sigma_ss, p_subset, T, pl)

nsubsets = size(sigma_ss,1);
rho = sigma_ss/sigma_0;
%gamma = (pl - T)./sqrt(1+rho.^2);
gamma = (pl - T)./sqrt(sigma_ss.^2+sigma_0^2);

pmd = zeros(nsubsets,1);
for i=1:nsubsets    
    pmd(i) = solSep.pl.findmaxF(rho(i),gamma(i));
end

phmi = p_subset'*pmd;


