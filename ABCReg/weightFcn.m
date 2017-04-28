function weights = weightFcn(distanceMaxt,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Calculate the weights using the Epanechnikov kernel
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weights = 1.5 * (epsilon)^(-1) * (1 - (distanceMaxt/epsilon).^(2));















