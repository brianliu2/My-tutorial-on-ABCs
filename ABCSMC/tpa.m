function t = tpa(epsilon,N,distance)
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Function is to decide the ration of samples remain
%   alive through the current iteration, while the living 
%   of sample depends on if the discrepancy
%   between the simulated dataset synthesized by sample
%   and real dataset is less than the current tolerance
%
%%%%%%%%%%%%%%%%%%%%%%%%%

%   Total number of samples
% N = size(x1,2);

% %   Distance between the synthetic dataset and the real dataset
% d = distanceFun_Update(x1,x2,x3,x4,x5,x6,y);

%   Vector contains indicators for implying if the samples
%   is alive
npa=sum((distance<epsilon),1); 

%   Scalar presents how many samples are alive during current iteration
pa=(npa>0);

%   Compute the ratio
t=sum(pa)/N;