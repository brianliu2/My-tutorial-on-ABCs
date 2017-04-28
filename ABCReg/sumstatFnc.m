function sumStat = sumstatFnc(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   function for calculating the summary statistics of pseudo-observations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



sumStat(1,:) = mean(x(1,:));
% sumStat(2,:) = var(x(1,:));
% sumStat(3,:) = var(x(1,:));
sumStat(2,:) = mean(x(2,:));
sumStat(3,:) = mean(x(3,:));
% sumStat(4,:) = var(x(1,:));
% sumStat(5,:) = var(x(2,:));
% sumStat(6,:) = var(x(3,:));

