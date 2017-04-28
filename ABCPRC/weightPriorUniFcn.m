function weight = weightPriorUniFcn(x,xHat,priorLow,priorUp,transSigma,weightPre)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This is the weight calculation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numerator = 1;
denominator = 0;

for i = 1:size(x,1)
    numerator = numerator * UniPDF(x(i,:),priorLow(i),priorUp(i));
    for n = 1:size(xHat,2)
        denominator = denominator + ...
                               weightPre(:,n)* normpdf(x(i,:),xHat(i,n),transSigma(i));
    end
end

weight = numerator / denominator;










