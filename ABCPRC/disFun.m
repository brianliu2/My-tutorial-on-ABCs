function distance = disFun(pseudoObs,Obs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   function for calculating the distance between the pseudo-observations
%   and the true observations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = 1:size(Obs,1)
    if max(pseudoObs(d,:)) > max(Obs(d,:))
        scalarFactor = max(pseudoObs(d,:));
        pseudoObs(d,:) = pseudoObs(d,:)/scalarFactor;
        Obs(d,:) = Obs(d,:)/scalarFactor;
    else
        scalarFactor = max(Obs(d,:));
        pseudoObs(d,:) = pseudoObs(d,:)/scalarFactor;
        Obs(d,:) = Obs(d,:)/scalarFactor;
    end    
end

distance = sqrt(sum(sum((pseudoObs - Obs).^2,2),1));





























