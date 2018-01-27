function fitDeposit(fl)

ctr = [.001 .005 .01 .05 .1 .5 1 5 10 25 50 75 100 250 500 750 1000];

[X,Y,Z] = plotT2(fl);

C           = contourf(X,Y,10.^Z, ctr);
[Cx,Cy,Cz]  = C2xyz(C);
Ca          = zeros(size(Cz));

for i = 1:length(Cz)
    Ca(i) = polyarea(Cx{i}, Cy{i})/1e6;
end

tephraFits(sqrt(Ca), Cz,  'exponential', 'dataType', 'isomass', 'BIS', [3,7])