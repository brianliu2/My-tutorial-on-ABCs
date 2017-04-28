%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Description:    Implementation is for running ABC-SMC for HS Model, in
%                         this code, parameter K_{d} and alpha_{d} are
%                         unknown
%
%   Input:
%                        paraSamples:     Initial samples for unknown parameters
%                        M:                       Integer factor
%                        alphaSMC:          Factor for weighting scheme
%                        N:                        Number of samples used
%                        delta:                   Discount factor for the transition
%                                                    kernel.
%                       options:                Options defined in the main
%                                                    function for running
%                                                    the ABC-SMC
%
%   Output:         
%                       paraSpace:  3-D matrix M*N*T contains estimations
%                                           at every iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
echo off;
clc

M                  = 3;
numParticles = 3;
alphaSMC     = 0.99;
epsilon          = 500;
epsilonTarget = 10;
essthresh = numParticles/2;

%   Number of processors required for parallel computing
Num_processor = 4;

%   Activate the parallel session
matlabpool('open', Num_processor);

%   Transition settings of parameter
delta              = 0.99;
a                    = (3*delta - 1)/(2*delta);
covParaCoeff  = 1 - a^2;

%   Regime definition for generating dataset
stateDim = 3;
paraDim  = 3;
timeLength = 100;
deltaT = 0.2;
DataSamples = timeLength/deltaT;
tspan = linspace(0,timeLength,DataSamples);

%   Definition of state and parameter settings
kd = 3; ad = 0.015; a0 = 0.03;
paraTrue = [kd,ad,a0];
xInit = [0 0 0];

%   Pre-allocate the memory for storing dataset
xSMC = zeros(M,stateDim,DataSamples);
weights = ones(1,numParticles)./numParticles;
epsilon_hist(1) = epsilon;

%   Generate particles for parameters at the first time step
lowBound = [2.9,0.01,0.025];
upBound  = [3.1,0.02,0.035];
paraSMC = zeros(paraDim,numParticles);
for d = 1:paraDim
    paraSMC(d,:) = lowBound(d) + (upBound(d) - lowBound(d)) * rand(1,numParticles);
end

paraSpace(:,:,1) = paraSMC;

%%  Generate synthetic data

options   = odeset('RelTol', 1e-02);
sol = ode45(@(t,y) hs_ode_KdAdA0(t,y,paraTrue),[0 timeLength],xInit,options);
xTrue = deval(sol,tspan);


%%  Main loop

%   Simulate the synthetic dataset by using the initial parameter particles
Loop = 1;
counter = 0;
xParallelLoop = generateData(M,numParticles,timeLength,tspan,xInit,paraSMC,options,Loop,epsilon);

while (epsilon > epsilonTarget)
    %   Compute the distance function
    distance = distanceFun(xParallelLoop(:,:,1,:),xParallelLoop(:,:,2,:),xParallelLoop(:,:,3,:),xTrue);
    
    % Calculate the next tolerance by following the condition presented in
    % the pseudo-code as PA() < alpha * PA()
    epsilon_old = epsilon;
    reflevel = alphaSMC * tpa(epsilon_old, numParticles, distance);
    epsilon  = fzero(@(epsilon) tpa(epsilon, numParticles, distance) - reflevel, 0, epsilon_old);
    epsilon_hist(Loop+1) = epsilon;       
    
    %   Determine if the current tolerance is less than the ideal tolerance,
    %   terminate the program if so
    if (epsilon < epsilonTarget)
        epsilon = epsilonTarget;
        epsilon_hist(Loop+1) = epsilon;
    end        
    
    %   Compute alive indicator vector considering the previous tolerance
    npa_old = sum((distance < epsilon_old),1);
    
    %   Compute alive indicator vector considering the current tolerance
    npa = sum((distance < epsilon),1);
    
    %   Find the indices for the alive or inactive samples 
    actIndex = find(npa_old > 0);
    unactIndex = find(npa_old == 0);
    
    %   Weight calculation according to the equations stated in the
    %   pseudo-code of ABC-SMC
    weights(1,actIndex) = weights(1,actIndex) .* npa(1,actIndex)./npa_old(1,actIndex);
    
    %   Assign the zero weights to inactive samples
    weights(1,unactIndex) = zeros(1,length(unactIndex));
    
    %   Normalized the weights
    weights = weights./sum(weights);     
        
    %   Resample if neccessary
    if (sum(weights.^2) * essthresh > 1)
        N_sons = rsdet(weights);
         [paraSMC,xParallelLoop] = copy(paraSMC,xParallelLoop,N_sons);       
        counter = counter + 1;
        counter
    end   
    
    %   Move the samples according to the transition kernel which is
    %   presented as Eqn.8 in the manuscript
    for c = 1:size(paraSMC,1)
            paraSMC(c,:) = a * paraSMC(c,:)...
                                    + ones(1,numParticles)*(1-a) * mean((paraSpace(c,:,Loop)'))...
                                    + covParaCoeff*std((paraSpace(c,:,Loop)'))...
                                    * randn(1,numParticles);
    end    
    
    %   Generate the synthetic dataset by using updated samples for
    %   parameters
    xParallelLoop = generateData(M,numParticles,timeLength,tspan,xInit,paraSMC,options,Loop,epsilon);
    
    Loop = Loop + 1;   
    paraSpace(:,:,Loop) = paraSMC;
    save('ABCSMC_HS.mat');
end

matlabpool('close')
Number_Workers = matlabpool('size');


paraSMC = paraSpace(:,:,end);
figure, 
plot(paraSMC');

figure,
subplot(3,1,1)
hist(paraSpace(1,:,end),40)
ylabel('density','FontSize',20);
xlabel('k_{d}, true: 3','FontSize',20)
set(gca,'FontSize',20);
title('estimate of K_{d}','FontSize',20);

subplot(3,1,2)
hist(paraSpace(2,:,end),40)
ylabel('density','FontSize',20);
xlabel('\alpha_{d}, true: 0.015','FontSize',20)
set(gca,'FontSize',20);
title('estimate of \alpha_{d}','FontSize',20);

subplot(3,1,3)
hist(paraSpace(3,:,end),40)
ylabel('density','FontSize',20);
xlabel('\alpha_{0}, true: 0.03','FontSize',20)
set(gca,'FontSize',20);
title('estimate of \alpha_{d}','FontSize',20);




































