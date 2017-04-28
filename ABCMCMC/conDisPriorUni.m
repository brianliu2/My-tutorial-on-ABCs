function pdf = conDisPriorUni(x,lowBound,upBound)

pdf = 1;

for i = 1:size(x,1)
    pdf = pdf * UniPDF(x(i,:),lowBound(i),upBound(i));
end