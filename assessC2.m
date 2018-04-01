% Second attempt with smaller variation

% clear
% 
% load C2
% 
% C = 100:100:900;
% 
% nhgt = length(data.conf.plumeHeight);
% ndur  = length(data.conf.duration);
% 
% stor = cell(nhgt,ndur);
% mT2 = zeros(numel(data.runlist.plumeHeight),1);
% 
% for i = 1:numel(data.runlist.plumeHeight)
%     disp(i)
%     
%     ihgt    = data.conf.plumeHeight == data.runlist.plumeHeight(i);
%     idur    = data.conf.duration == data.runlist.duration(i);
%     
%     [massT2, massfit] = fitDeposit(['C2/OUT/', num2str(i),'.out'], C);
%     if ~isnan(massT2)
%         stor{ihgt, idur}  = massfit;
%         mT2(i) = massT2;
%     end
% end
% 
% storTot = zeros(nhgt,ndur,length(C),3);
% count = 1;
% leg = cell(length(C),1);
% for iH = 1:nhgt
%     for iD = 1:ndur
%         for iC = 1:length(C)
%             storT(iH,iD,iC,1) = stor{iH,iD}{iC}.powerlaw.mass_kg;
%             storT(iH,iD,iC,2) = mT2(count);
%             storT(iH,iD,iC,3) = storT(iH,iD,iC,1) / storT(iH,iD,iC,2);
%             leg{iC} = num2str(C(iC));
%         end
%         count = count+1;
%     end
% end

figure;
% a(1) = subplot(2,1,1); hold on, box on
% a(2) = subplot(2,1,2); hold on, box on

for iD = 1:ndur
    subplot(2,1,iD); hold on, box on
    for iH = 1:nhgt
        plot(C, squeeze(storT(iH,iD,:,3)))
        
        
    end
end
legend(sprintfc('%d',data.conf.plumeHeight./1e3));

