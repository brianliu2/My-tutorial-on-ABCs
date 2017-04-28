function stateData = generateData(M,numParticles,timeLength,tspan,xInit,paraSMC,options,loop,epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Generate synthetic dataset by using the current parameter particles 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:M
    parfor j = 1:numParticles
        %   solve the ODEs of LV system
        solSMC = ode45(@(t,y) hs_ode_KdAdA0(t,y, paraSMC(:,j)),[0 timeLength],xInit,options);        
        %   Use the last column in the solution as the underlying
        %   state, if the solver is underflow.
        if(solSMC.x(end) >= timeLength)
            temp = deval(solSMC,tspan);
        else
            temp = deval(solSMC,tspan);
        end

        stateData(i,j,:,:) = temp;
        fprintf('%d particles in %d M of %d loop and recent epsilon is %d \n',j,i,loop,epsilon);        
    end
end