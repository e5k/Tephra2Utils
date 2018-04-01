load C4
fl = dir('C3/OUT/*.out');
storA = zeros(length(fl),2);
storM = zeros(length(fl),2);
param = struct;

% for i = 1:length(fl)
%     display(i)
%     [~,~,~,A] = plotT2(['C3/OUT/', num2str(i), '.out'], '-noplot');
%     storA(i,:) = [A(2,1),A(3,1)];
%     storM(i,:) = [A(2,2),A(3,2)];
% end
% save('storAC3.mat','storA','storM')

load storAC3
%% check data.master
figure
hold on

ax1 = subplot(2,2,1); hold on, box on
xlabel('Magnitude');
ylabel('1 kg/m^2 isomass area^2 (km)');
xlim([1.5 6.6]);

ax2 = subplot(2,2,2); hold on, box on
xlabel('Magnitude');
ylabel('10 kg/m^2 isomass area^2 (km)');
xlim([1.5 6.6]);

ax3 = subplot(2,2,3); hold on, box on
xlabel('Magnitude');
ylabel('% of total mass');
xlim([1.5 6.6]);

ax4 = subplot(2,2,4); hold on, box on
xlabel('Magnitude');
ylabel('% of total mass');
xlim([1.5 6.6]);

dur = [1,6,12];
med = [-1,1];

durS = lines(3);
medS = {'^','o'};
for iD = 1:length(data.conf.duration)
    for iM = 1:length(data.tgsd.med)
        idx = data.master.duration == dur(iD) & data.master.med == med(iM);
        axes(ax1);
        if iM == 1
            plot(log10(data.master.mass(idx))-7, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(log10(data.master.mass(idx))-7, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:),'MarkerFaceColor', durS(iD,:))
        end
        
        axes(ax2);
        if iM == 1
            plot(log10(data.master.mass(idx))-7, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(log10(data.master.mass(idx))-7, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', durS(iD,:))
        end
            
        axes(ax3);
        if iM == 1
            plot(log10(data.master.mass(idx))-7, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(log10(data.master.mass(idx))-7, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
        end
            
        axes(ax4);
        if iM == 1
            plot(log10(data.master.mass(idx))-7, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(log10(data.master.mass(idx))-7, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
        end
    end
end

ylim(ax3,[60 100]);
ylim(ax4,[60 100]);















figure
hold on

ax1 = subplot(2,2,1); hold on, box on
xlabel('Plume height (km)');
ylabel('1 kg/m^2 isomass area^2 (km)');

ax2 = subplot(2,2,2); hold on, box on
xlabel('Plume height (km)');
ylabel('10 kg/m^2 isomass area^2 (km)');

ax3 = subplot(2,2,3); hold on, box on
xlabel('Plume height (km)');
ylabel('% of total mass');

ax4 = subplot(2,2,4); hold on, box on
xlabel('Plume height (km)');
ylabel('% of total mass');

dur = [1,6,12];
med = [-1,1];

durS = lines(3);
medS = {'^','o'};
for iD = 1:length(data.conf.duration)
    for iM = 1:length(data.tgsd.med)
        idx = data.master.duration == dur(iD) & data.master.med == med(iM);
        axes(ax1);
        if iM == 1
            plot(data.master.plumeHeight(idx)./1e3, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(data.master.plumeHeight(idx)./1e3, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:),'MarkerFaceColor', durS(iD,:))
        end
        
        axes(ax2);
        if iM == 1
            plot(data.master.plumeHeight(idx)./1e3, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(data.master.plumeHeight(idx)./1e3, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', durS(iD,:))
        end
            
        axes(ax3);
        if iM == 1
            plot(data.master.plumeHeight(idx)./1e3, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(data.master.plumeHeight(idx)./1e3, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
        end
            
        axes(ax4);
        if iM == 1
            plot(data.master.plumeHeight(idx)./1e3, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
        else
            plot(data.master.plumeHeight(idx)./1e3, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
        end
    end
end

ylim(ax3,[60 100]);
ylim(ax4,[60 100]);


% 
% figure
% hold on
% 
% ax1 = subplot(2,2,1); hold on, box on
% xlabel('Magnitude');
% ylabel('1 kg/m^2 isomass area^2 (km)');
% 
% ax2 = subplot(2,2,2); hold on, box on
% xlabel('Magnitude');
% ylabel('10 kg/m^2 isomass area^2 (km)');
% 
% ax3 = subplot(2,2,3); hold on, box on
% xlabel('Magnitude');
% ylabel('% of total mass');
% 
% ax4 = subplot(2,2,4); hold on, box on
% xlabel('Magnitude');
% ylabel('% of total mass');
% 
% dur = [1,6,12];
% med = [-1,1];
% 
% durS = lines(3);
% medS = {'^','o'};
% for iD = 1:length(data.conf.duration)
%     for iM = 1:length(data.tgsd.med)
%         idx = data.master.duration == dur(iD) & data.master.med == med(iM);
%         axes(ax1);
%         if iM == 1
%             plot(log10(data.master.MER(idx))+3, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
%         else
%             plot(log10(data.master.MER(idx))+3, sqrt(storA(idx,1)), 'o', 'Color', durS(iD,:),'MarkerFaceColor', durS(iD,:))
%         end
%         
%         axes(ax2);
%         if iM == 1
%             plot(log10(data.master.MER(idx))+3, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
%         else
%             plot(log10(data.master.MER(idx))+3, sqrt(storA(idx,2)), 'o', 'Color', durS(iD,:), 'MarkerFaceColor', durS(iD,:))
%         end
%             
%         axes(ax3);
%         if iM == 1
%             plot(log10(data.master.MER(idx))+3, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
%         else
%             plot(log10(data.master.MER(idx))+3, storM(idx,1)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
%         end
%             
%         axes(ax4);
%         if iM == 1
%             plot(log10(data.master.MER(idx))+3, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor', 'w')
%         else
%             plot(log10(data.master.MER(idx))+3, storM(idx,2)./data.master.mass(idx)'.*100, 'o', 'Color', durS(iD,:), 'MarkerFaceColor',  durS(iD,:))
%         end
%     end
% end