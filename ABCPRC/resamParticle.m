function reParticles = resamParticle(weights,preParticles)

reParticles = zeros(size(preParticles,1), size(preParticles,2));
for n = 1:size(preParticles,2)
    for r = 1:size(preParticles,1)
        u = rand;
        qsum = 0;
        reFlag = 1;
        indx = 1;
        while indx < size(preParticles,2) && reFlag == 1
            qsum = qsum + weights(indx);
            if qsum > u
                reParticles(r,n) = preParticles(r,indx);
                reFlag = 0;
            end
            indx = indx + 1;
        end
    end
end