load C4
fl = dir('C3/OUT/*.out');

stor = struct;
cnt = [1,5,10,25,50,75,100,150,200,300,400,500,750,1000,1500,2000,3000,4000,5000,7500,10000];
for i = 1:length(fl)
    display(i)
    [~,~,~,stor(i).area,stor(i).mass,stor(i).cnt] = plotT2(['C3/OUT/', num2str(i), '.out'], 'contours', cnt, '-noplot');
    stor(i).massT = data.master.mass(i);
    stor(i).fit   = tephraFits(sqrt(stor(i).area), stor(i).cnt, 'exponential','segments', [1 4],'dataType','isomass');
    stor(i).mass  = stor(i).fit.exponential.mass_kg;
    print(gcf, ['C3/FIG/', num2str(i),'.pdf'], '-dpdf');
    close(gcf);
end
save('fit_exp.mat', 'stor');

figure,
plot([0,max([stor.mass])], [0,max([stor.mass])], '-k'); hold on
plot([stor.massT],[stor.mass], '.r', 'markersize',12)
xlabel('Tephra2 mass');
ylabel('Fit mass');
axis square equal
axis([0,max([stor.mass]),0,max([stor.mass])])