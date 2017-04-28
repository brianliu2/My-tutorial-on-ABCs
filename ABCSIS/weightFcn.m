function weight = weightFcn(x,priorMu,priorSigma,transMu,transSigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is the weight calculation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numerator = 1;
denominator = 1;

for d = 1:size(x,1)
    numerator = numerator * normpdf(x(d),priorMu(d),priorSigma(d))* normpdf(transMu(d),x(d),transSigma(d));
    denominator = denominator * normpdf(transMu(d),priorMu(d),priorSigma(d)) * normpdf(x(d),transMu(d),transSigma(d));
end

weight = numerator / denominator;