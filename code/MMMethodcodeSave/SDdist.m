function [ zscorecutoff,EV,Var,SD ] = SDdist( x ,confidencebound)
%x is a vector of numbers
%this computes the expected value of x, the variance of x, the 
%distribution by standard deviations, and the confidence interval. 
%
EV = sum(x)/length(x);
Var = sum((x-EV).^2)/length(x);
SD = sqrt(Var);
zscore = (x-EV)/SD;

sortedabszscore = sort(abs(zscore));
cutoff = ceil(confidencebound*length(x));
zscorecutoff = sortedabszscore(cutoff);





end

