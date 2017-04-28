%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   ABC-MCMC on LV model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
echo off;
% close all;
clc

global lowerBound;
global upperBound;
global priorMu;
global priorSigma;
global sigmaTrans;


%   Definition of algorithm
MHIteration = 2000;
epsilon = 10;
acceptCont = 0;


%   Initialization of generating data
timeLength = 100;
deltaT = 0.2;
dataPoint = timeLength/deltaT;
tspan = linspace(0,timeLength,dataPoint);
xInit = [0 0 0];
options = odeset('RelTol',1e-02);

%   Initialization of particles
lowerBound = [0;0];
upperBound = [10;0.1];
paraDim = 2;
paraAccept = zeros(paraDim,MHIteration);
% priorMu = 0.7 * ones(2,1);
% priorSigma = 0.01*ones(2,1);
% sigmaTrans = 0.023 * ones(paraDim,1);
% sigmaTrans = [0.01;0.001];
sigmaTrans = [0.0001;0.000007];
%%  Generate true observations
kd = 3; ad = 0.015;
paraTrue = [kd,ad];
sol = ode45(@(t,y) hs_ode_KdAd(t,y,paraTrue),[0 timeLength],xInit);
xTrue = deval(sol,tspan);
% xTrue = xTrue + sqrt(0.0005) * randn(size(xTrue));

%%  Main Loop
Loop = 2;
startMCMCTime = cputime;
distMatrix  = zeros(1,MHIteration);

%   pre-allocate the memory for storing particles
paraMHSpace = zeros(paraDim,MHIteration);

%   generate the initial particles for parameters
paraMH = [3.5;0.02];
% paraMH = lowerBound + (upperBound - lowerBound) .* rand(paraDim,1);
% paraMH = priorMu + sqrt(priorSigma).* randn(paraDim,1);
paraMHSpace(:,1) = paraMH;

while (Loop-1 < MHIteration)
    %   1. perturb the current particles 
    paraMH_hat = paraMH + sqrt(sigmaTrans).* randn(size(paraMH));
    
    %   2. synthesize pseudo-observation
    solMH = ode45(@(t,y) hs_ode_KdAd(t,y,paraMH_hat),[0 timeLength],xInit);
    xMH = deval(solMH,tspan);
%     xMH = xMH + sqrt(0.0005) * randn(size(xMH));
    
    %   3. calculate the distance
    distance = disFun(xMH,xTrue);
    
    %   4. accept the current particles if distance is less than the
    %       tolerance 
    if distance < epsilon
        acceptCont = acceptCont +1;
%         if acceptCont == 100
%             Loopfor100Particles = Loop;
%         end
        TarDis1 = conDisPriorUni(paraMH_hat,lowerBound,upperBound)+conDis(paraMH_hat,paraMH,sigmaTrans);
        TarDis2 = conDisPriorUni(paraMH,lowerBound,upperBound)+conDis(paraMH,paraMH_hat,sigmaTrans);
        if log(rand) < (TarDis1 - TarDis2)
            paraMH = paraMH_hat;
            paraMHSpace(:,Loop) = paraMH;
        else
            paraMHSpace(:,Loop) = paraMH;
        end
        disp([num2str(acceptCont),' particles are accpeted within ', num2str(Loop-1), 'MH loops']);
    else
        paraMHSpace(:,Loop) = paraMH;
    end
    distMatrix(:,Loop) = distance;
    Loop = Loop + 1;
    Loop
end


acceptRate = acceptCont/MHIteration;
finishMCMCTime = cputime - startMCMCTime;

%%  Plot

% [fpara1, xpara1] = hist(paraAccept(1,:),numParticles);
% para1PDF = normpdf(paraAccept(1,:),mean(paraAccept(1,:)),var(paraAccept(1,:)));

% % subplot(1,2,2)
% % [fpara2, xpara2] = hist(paraAccept(2,:),40);

%   contour

% mu = mean(paraAccept')';
% sigma = cov(paraAccept');

% a1 = min(paraAccept(1,:)); b1 = max(paraAccept(1,:)); 
% a2 = min(paraAccept(2,:)); b2 = max(paraAccept(2,:));
% mvg(mu,sigma,a1,b1,a2,b2,100);

paraAcceptMH = paraMHSpace(:,2000:end);

muMH = mean(paraAcceptMH')';
sigmaMH = cov(paraAcceptMH');

a1MH = min(paraAcceptMH(1,:)); b1MH = max(paraAcceptMH(1,:)); 
a2MH = min(paraAcceptMH(2,:)); b2MH = max(paraAcceptMH(2,:));
mvg(muMH,sigmaMH,a1MH,b1MH,a2MH,b2MH,100);
plot(0.5,0.5,'+');
%   3-D plot
% figure
% bar(xpara1,fpara1/trapz(xpara1,fpara1));
% hold on;
% plot(xpara1,para1PDF,'r');























