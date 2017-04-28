function regCoef = regCoefFnc(sumStat,sumStatTrue,weights,paraAccept)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculate the regression coefficents matrix 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Replicate the summary statistics of true observations in order to meet
%   the dimension requirment of summary statistics of pseudo-observation
sumStatTrue = sumStatTrue(:,ones(1,size(sumStat,2)));

%   Differences of summary statistics between pseudo-observations and
%   observations
sumStatDif = sumStat - sumStatTrue;

%   Calculate the regression coefficients matrix
regCoef = inv(sumStatDif * diag(weights) * sumStatDif') * sumStatDif * diag(weights) * paraAccept';













