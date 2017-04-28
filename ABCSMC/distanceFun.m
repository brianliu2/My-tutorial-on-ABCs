function d = distanceFun(x1,x2,x3,trueObs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The distance function for measuring distance the 
%   pseudo observations and true observations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Pre-allocate the memory for storing distances
%   the dimension of distance matrix is M * N, in which
%   M is the integer factor and N is the number of particles

d = zeros(size(x1,1),size(x1,2));

for i = 1:size(x1,2)    %   Loop for number of particles N    
    for j = 1:size(x1,1)    %   Loop for integer factor M
        %   1. reshape the dataset
        pseudoObs1 = reshape(x1(j,i,:,:),1,size(trueObs,2));
        pseudoObs2 = reshape(x2(j,i,:,:),1,size(trueObs,2));
        pseudoObs3 = reshape(x3(j,i,:,:),1,size(trueObs,2));
        trueObs1 = trueObs(1,:);
        trueObs2 = trueObs(2,:);
        trueObs3 = trueObs(3,:);
        
        %   2. scale down the dataset
        if max(pseudoObs1) > max(trueObs1)
            scaleFactor = max(pseudoObs1);
            pseudoObs1 = pseudoObs1/scaleFactor;
            trueObs1 = trueObs1/scaleFactor;
        else
            scaleFactor = max(trueObs1);
            pseudoObs1 = pseudoObs1/scaleFactor;
            trueObs1 = trueObs1/scaleFactor;
        end
        
        if max(pseudoObs2) > max(trueObs2)
            scaleFactor = max(pseudoObs2);
            pseudoObs2 = pseudoObs2/scaleFactor;
            trueObs2 = trueObs2/scaleFactor;
        else
            scaleFactor = max(trueObs2);
            pseudoObs2 = pseudoObs2/scaleFactor;
            trueObs2 = trueObs2/scaleFactor;
        end
        
        if max(pseudoObs3) > max(trueObs3)
            scaleFactor = max(pseudoObs3);
            pseudoObs3 = pseudoObs3/scaleFactor;
            trueObs3 = trueObs3/scaleFactor;
        else
            scaleFactor = max(trueObs3);
            pseudoObs3 = pseudoObs3/scaleFactor;
            trueObs3 = trueObs3/scaleFactor;
        end
        
        %   3. the dimension of true state is 2-D
        %       y(state dimension, data points for each time interval).
        %       the dimension of synthetic data is 4-D
        %       x(integer factor, number of particles, state dimension, data points for each time interval)
        %       In calculation, the first summation is for summing
        %       all states in the system, the second summation is 
        %       for data points each state in each time interval.
        d(j,i) = sqrt(sum(sum((pseudoObs1 - trueObs1).^2) + sum((pseudoObs2 - trueObs2).^2)+...
                          sum((pseudoObs3 - trueObs3).^2)));
    end
end








