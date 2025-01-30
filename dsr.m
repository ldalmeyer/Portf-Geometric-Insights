function [deflatedSharpe, pvalue] = dsr(N, T, observedsharpe, variancesharpe, skew, kurtosis)
%nmatare quanttools (github)

   expectedsharpe=0; %the null hypothesis that strategies are not better than 0
   gamma = -psi(1); %E-M constant
   %compute the expected maximum given a sample size of N and i.i.d assumptions
   maxn = (1 - gamma) * norminv(1 - (1 / N)) + gamma * norminv(1 - (1 / N) * exp(-1)); %SR_0
   nullsharpe = (expectedsharpe + sqrt(variancesharpe) * maxn);
   
   numerator = (observedsharpe - nullsharpe) * sqrt(T - 1);
   adjskew = 1 - skew * observedsharpe';
   adjkurt = ((kurtosis - 1) / 4) * observedsharpe' .^ 2;
    
    % this means that (adj.skew + adj.kurt) < 0, and sqrt(# < 0) == Nan
    if (adjskew + adjkurt) < 0 
        denominator = 1 / 1e6; % will force 1 - pvalue to 1L
    else 
        denominator = sqrt(adjskew + adjkurt);
    end
     
    pvalue = (1 - normcdf(numerator / denominator));
    deflatedSharpe =  (normcdf(numerator / denominator) .* observedsharpe) + (pvalue .* expectedsharpe);   
   
end