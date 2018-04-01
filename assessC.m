clear

load C_tmp

C = 100:100:800;

nwind = length(data.wind.vel);
nhgt = length(data.conf.plumeHeight);
ndur  = length(data.conf.duration);
ntgsd = length(data.tgsd.med);

stor = cell(nhgt,nwind,ndur,ntgsd);

for i = 257:numel(data.runlist.plumeHeight)
    disp(i)
    
    ihgt    = data.conf.plumeHeight == data.runlist.plumeHeight(i);
    iwind   = data.wind.vel == data.runlist.vel(i);
    idur    = data.conf.duration == data.runlist.duration(i);
    itgsd   = data.tgsd.med == data.runlist.med(i);
    
    [massT2, massfit] = fitDeposit(['C/OUT/', num2str(i),'.out'], C);
    if ~isnan(massT2)
        stor{ihgt, iwind, idur, itgsd}  = massfit.powerlaw.mass_kg;
    end
end
