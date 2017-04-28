function prob = conDis(x, mean,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   Function for calculating the acceptance ratio in the
%   Metropolis-Hastings algorithm
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prob = 0;

for d = 1:size(x,1)
    prob = prob + -0.5 * (x(d)-mean(d))^2/sigma(d);
end

% % % global lowerBound;
% % % global upperBound;
% % % global priorMu;
% % % global priorSigma;
% % % global sigmaTrans;
% % % numerator = pdfFN(paraMH_hat,lowerBound,upperBound,paraMH,sigmaTrans,'unif')...
% % %                  * pdfFN(paraMH_hat,lowerBound,upperBound,paraMH,sigmaTrans,'normal');
% % % denominator = pdfFN(paraMH,lowerBound,upperBound,paraMH_hat,sigmaTrans,'unif')...
% % %                  * pdfFN(paraMH,lowerBound,upperBound,paraMH_hat,sigmaTrans,'normal');
% % % prob = numerator / denominator;
% % 
% % % numerator = pdfFN(paraMH_hat,lowerBound,upperBound,priorMu,priorSigma,'normal')...
% % %                  * pdfFN(paraMH_hat,lowerBound,upperBound,paraMH,sigmaTrans,'normal');
% % % denominator = pdfFN(paraMH,lowerBound,upperBound,priorMu,priorSigma,'normal')...
% % %                  * pdfFN(paraMH,lowerBound,upperBound,paraMH_hat,sigmaTrans,'normal');
% % % prob = numerator / denominator;
% % numerator = log(pdfFN(paraMH_hat,lowerBound,upperBound,priorMu,priorSigma,'normal'))...
% %                  + log(pdfFN(paraMH_hat,lowerBound,upperBound,paraMH,sigmaTrans,'normal'));
% % denominator = log(pdfFN(paraMH,lowerBound,upperBound,priorMu,priorSigma,'normal'))...
% %                  + log(pdfFN(paraMH,lowerBound,upperBound,paraMH_hat,sigmaTrans,'normal'));
% % prob = numerator - denominator;
% % 
% % 
% % % numratorPrior = zeros(size(paraMH,1),1);
% % % numratorTran = zeros(size(paraMH,1),1);
% % % denominatorPrior = zeros(size(paraMH,1),1);
% % % denominatorTran = zeros(size(paraMH,1),1);
% % % 
% % % prob = 1;
% % % 
% % % for d = 1:size(paraMH_hat)
% % %     prob
% % % end
% % 
% % 









