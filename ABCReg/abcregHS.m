clear all
clc

%   Definition of algorithm
numParticles = 1000;
epsilon = 0.015;
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
xTrueScaled = scaleFcn(xTrue);
sumStatTrue = sumstatFnc(xTrueScaled);

%%  Main Loop
Loop = 0;
startRegTime = cputime;
% sumStat = zeros(6,numParticles);
sumStat = zeros(3,numParticles);
distMatrix  = zeros(1,numParticles);

while (acceptCont < numParticles)
    %   1. generate particle for parameter
    paraRej = lowBound + (upBound - lowBound) .* rand(paraDim,1);
    
    %   2. synthesize pseudo-observation
    para = [paraRej(1);paraRej(2)];
    solRej = ode45(@(t,y) hs_ode_KdAd(t,y,para),[0 timeLength],xInit);
    xRej = deval(solRej,tspan);
    pseudoX = scaleFcn(xRej);
    
    %   3. calculate the distance
    pseudosumStat = sumstatFnc(pseudoX);
    distance = disFun(pseudosumStat,sumStatTrue);
    
    %   4. accept the current particles if distance is less than the
    %       tolerance 
    if distance < epsilon
        acceptCont = acceptCont +1;
        paraAccept(:,acceptCont) = paraRej;
        sumStat(:,acceptCont) = pseudosumStat;
        distMatrix(:,acceptCont) = distance;
        disp([num2str(acceptCont),' particles are accpeted within ', num2str(Loop), 'loops']);
        save('ABCReg_HS.mat');
    end    
    
    Loop = Loop + 1;
    disp([num2str(acceptCont),' particles are accpeted within ', num2str(Loop), ' loops']);    
end

save('ABCReg_HS.mat');
%   5. Update the accepted particles according to the local regerssion
%   adjustment
weights = weightFcn(distMatrix,epsilon);
regCoef = regCoefFnc(sumStat,sumStatTrue,weights,paraAccept);

%   Replicate the summary statistics of the true observations
sumStatTrueRep = sumStatTrue(:,ones(1,size(sumStat,2)));

%   Update the accepted particles
paraAcceptReg = paraAccept - regCoef' * (sumStat-sumStatTrueRep);

acceptRate = numParticles/Loop;
finishRegTime = cputime - startRegTime;

mu = mean(paraAccept')';
sigma = cov(paraAccept');

a1 = min(paraAccept(1,:)); b1 = max(paraAccept(1,:)); 
a2 = min(paraAccept(2,:)); b2 = max(paraAccept(2,:));
mvg(mu,sigma,a1,b1,a2,b2,100);
plot(3,0.015,'+');

muReg = mean(paraAcceptReg')';
sigmaReg = cov(paraAcceptReg');

a1Reg = min(paraAcceptReg(1,:)); b1Reg = max(paraAcceptReg(1,:)); 
a2Reg = min(paraAcceptReg(2,:)); b2Reg = max(paraAcceptReg(2,:));
mvg(muReg,sigmaReg,a1Reg,b1Reg,a2Reg,b2Reg,100);
plot(3,0.015,'+');











































