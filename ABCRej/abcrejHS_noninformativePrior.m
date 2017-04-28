%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   ABC-rejection on LV model
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc

%   Definition of algorithm
numParticles = 1000;
epsilon = 0.7;
acceptCont = 0;

%   Initialization of generating data
timeLength = 100;
deltaT = 0.2;
dataPoint = timeLength/deltaT;
tspan = linspace(0,timeLength,dataPoint);
xInit = [0 0 0];
options = odeset('RelTol',1e-02);

%   Initialization of particles
lowBound = [0;0];
upBound  = [10;1];
paraDim = 2;
paraAccept = zeros(paraDim,numParticles);

%%  Generate true observations
kd = 3; ad = 0.015;
paraTrue = [kd;ad];
sol = ode45(@(t,y) hs_ode_KdAd(t,y,paraTrue),[0 timeLength],xInit);
xTrue = deval(sol,tspan);

%%  Main Loop
Loop = 0;
startRejTime = cputime;
% load('ABCRej_HS.mat');
while (acceptCont < numParticles)
    %   1. generate particle for parameter
    paraRej = lowBound + (upBound - lowBound) .* rand(paraDim,1);
%     paraRej = priorMu + sqrt(priorSigma) * randn(paraDim,1);
    
    %   2. synthesize pseudo-observation
    para = [paraRej(1);paraRej(2)];
    solRej = ode45(@(t,y) hs_ode_KdAd(t,y,para),[0 timeLength],xInit);
    xRej = deval(solRej,tspan);
    
    %   3. calculate the distance
    distance = disFun(xRej,xTrue);
    
    %   4. accept the current particles if distance is less than the
    %       tolerance 
    if distance < epsilon
        acceptCont = acceptCont +1;
        paraAccept(:,acceptCont) = paraRej;
        disp([num2str(acceptCont),' particles are accpeted within ', num2str(Loop), 'loops']);
        save('ABCRej_HS2.mat');
    end
    
    Loop = Loop + 1;
    disp([num2str(acceptCont),' particles are accpeted within ', num2str(Loop), ' loops']);    
end
acceptRate = numParticles/Loop;
finishRejTime = cputime - startRejTime;
save('ABCRej_HS2.mat');

%%  Plot
%   histogram
% % figure,
% % subplot(1,2,1)
[fpara1, xpara1] = hist(paraAccept(1,:),numParticles);
para1PDF = normpdf(paraAccept(1,:),mean(paraAccept(1,:)),var(paraAccept(1,:)));

% % subplot(1,2,2)
% % [fpara2, xpara2] = hist(paraAccept(2,:),40);

%   contour
mu = mean(paraAccept')';
sigma = cov(paraAccept');

a1 = min(paraAccept(1,:)); b1 = max(paraAccept(1,:)); 
a2 = min(paraAccept(2,:)); b2 = max(paraAccept(2,:));
mvg(mu,sigma,a1,b1,a2,b2,100);


%   3-D plot
% figure
% bar(xpara1,fpara1/trapz(xpara1,fpara1));
% hold on;
% plot(xpara1,para1PDF,'r');























