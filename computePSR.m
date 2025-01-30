function [psr] = computePSR(data, sr_ref)

%mean_d = mean(data);
%sd_d = std(data);
skew_d = skewness(data);
kurt_d = kurtosis(data);
obs = size(data,1);

sr = sharpe_L(data);
psrStat_Num = ((sr - sr_ref) * (obs - 1) ^ 0.5);
psrStat_Den = (1 - (skew_d * sr) + ((kurt_d-1)/4 * sr^2)) ^ 0.5;
psr = normcdf(psrStat_Num/psrStat_Den);
