function distance = disFun(pseudoObs,Obs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   function for calculating the distance between the pseudo-observations
%   and the true observations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

distance = sqrt(sum(sum((pseudoObs - Obs).^2,2),1));





























