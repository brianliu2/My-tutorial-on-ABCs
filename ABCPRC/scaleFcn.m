function valScaled = scaleFcn(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   function for calculating the distance between the pseudo-observations
%   and the true observations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

valScaled = zeros(size(x));
for d = 1:size(x,1)
    x(d,:) = x(d,:)/max(x(d,:));
end