function [LL] = ARCHMneg_LL(param_vec,ddata)
%This function returns the negative effective log-likelihood given the parameters and
%data
mmu_s = param_vec(1);
rrho_s = param_vec(2);

% GARCH innovation parameters
oomega = param_vec(3)^2;
aalpha = param_vec(4);
ggamma = param_vec(5);

T = length(ddata);
innovs = zeros(1,T);
var_est = innovs;

innovs(1) = sqrt(oomega/(1-aalpha));
var_est(1) = innovs(1)^2;

for t=2:T
    var_est(t) = oomega + aalpha*innovs(t-1)^2;
    innovs(t) = ddata(t) - ( (1-rrho_s)*mmu_s + rrho_s*ddata(t-1) + ggamma*sqrt(var_est(t)) );
end

LL = .5*sum(log(var_est)+innovs.^2./var_est);

end

