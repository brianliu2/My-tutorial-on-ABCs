%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   ABC-MCMC on Heat shock model, this is the one for the most stiff
%   parameters k_{d} and alpha_{d}
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
echo off;
clc

randn('state',sum(100*clock))
%   Definitions for running ABC-PRC
numParticles = 2500;
SMCIteration = 239;
paraDim = 3;
stateDim = 3;
epsilon = linspace(488,10,239);
paraSMCSpace = zeros(paraDim,numParticles,SMCIteration);
%   Definitions for the prior distribution
lowBound = [2.2;0.01;0.02];
upperBound   = [3.2;0.02;0.04];
paraSMC = zeros(paraDim,numParticles);
for d = 1:paraDim
    paraSMC(d,:) = lowBound(d) + (upperBound(d) - lowBound(d)) * rand(1,numParticles);
end
paraSMCSpace(:,:,1) = paraSMC;
%%   Generate the true observations
timeLength = 100;
deltaT = 0.2;
samplePoints = timeLength/deltaT;
tspan = linspace(0,timeLength,samplePoints);
kd = 3; ad = 0.015; a0 = 0.03; as = 3; ks = 0.05; ku = 0.0254;
paraTrue = [kd,ad,a0,as,ks,ku];
xInit = [0,0,0];
options   = odeset('RelTol', 1e-02);
sol = ode45(@(t,y) hs_ode(t,y,paraTrue),[0 timeLength],xInit,options);
xTrue = deval(sol,tspan);


%%  ABC-SIS
xSMCSpace = zeros(stateDim,samplePoints,numParticles,SMCIteration);


%   Generate the pseudo-observations using the intial mu particles
xSMC                   = zeros(stateDim,samplePoints,numParticles);
xSMCSpace(:,:,:,1) = xSMC;
%   Define the weights for the initial particles
weights = ones(1,numParticles)/numParticles;
weightSpace = zeros(1,numParticles,SMCIteration);
weightSpace(:,:,1) = weights;
%   Define the sigma for the transition kernel
transSigma = [2;0.01;0.02];

resampleCnt = 0;
paraSMCHat = zeros(paraDim,numParticles);
evalTime = zeros(numParticles,SMCIteration);
%   Main Loop for ABCPRC
for t = 2:SMCIteration
    acceptCnt = 0;
    while acceptCnt ~= numParticles
        acceptFlag = 'false';
        evalCnt = 0;
        while (strcmp(acceptFlag,'false') == 1)
            %   1. select the particles associate the previous weights
            for d = 1:paraDim
                paraSMCHat(d,:) = randsample(paraSMCSpace(d,:,t-1),numParticles,true,weights);                
            end
            %   2. perturb the particles according to the transition kernel
            for d = 1:paraDim
                paraSMC(d,:) = paraSMCHat(d,:) + transSigma(d)* randn(1,numParticles);
            end
            
            %   3. generate the pseudo-observations by using a given particle
            paraCur = [paraSMC(1,acceptCnt+1),paraSMC(2,acceptCnt+1),paraSMC(3,acceptCnt+1),as,ks,ku];
            solSMC = ode45(@(t,y) hs_ode(t,y,paraCur),[0 timeLength],xInit,options);
            xSMC = deval(solSMC,tspan);
            %   4. Calculate the distance between x^{ast} and x
            distance = disFun(xSMC,xTrue);

            %   5. determine if the current particle should be accepted
            if (distance < epsilon(t-1))                
                weights(acceptCnt+1) = weightPriorUniFcn(paraSMC(:,acceptCnt+1),...
                                                  paraSMCSpace(:,:,t-1),lowBound,...
                                                  upperBound,transSigma,weightSpace(:,:,t-1))+eps;
                xSMCSpace(:,:,acceptCnt+1,t) = xSMC;
                paraSMCSpace(:,acceptCnt+1,t)=paraSMC(:,acceptCnt+1);
                acceptCnt = acceptCnt +1;
                acceptFlag = 'true';                
            end
            evalCnt = evalCnt + 1;
        end
        evalTime(acceptCnt,t) = evalCnt;
        disp([num2str(acceptCnt),' particles accepted in ',num2str(t), ' iteration']);
    end
    %   6. Normalize the weights
    weights = weights / sum(weights);
    weightSpace(:,:,t) = weights;
    %   6. resample the particles if a degeneracy problem is observed
%     if (sum(weights.^2)^(-1) < numParticles/2)
%         paraSMCSpace(:,:,t) = randsample(paraSMCSpace(:,:,t),numParticles,true,weights);
%         weights = ones(1,numParticles)/numParticles;
%         resampleCnt = resampleCnt + 1;
%     end
    save('ABCSIS_HS_ThreeStiff_informative.mat')
    t
end

% plotThesis(paraSMCSpace,numParticles,epsilon);







