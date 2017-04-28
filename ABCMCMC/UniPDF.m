function pdfUni = UniPDF(x,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Calculate the pdf of uniform distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%

if x>a && x<b
    pdfUni = 1/(b-a);
else
    pdfUni = 0;
end