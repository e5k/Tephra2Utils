
function TGSD = plotTGSD(mu, sigma)

range = -20:20;
dist  = normpdf(range, mu, sigma);
dist2 = normcdf(range, mu, sigma);

idx = dist>1e-3;

figure, 
bar(range(idx), dist(idx).*100);
xlabel('\phi')
ylabel('Wt. %')

yyaxis right
plot(range(idx), dist2(idx).*100, '-k', 'Linewidth',1);

TGSD = [range(idx)', dist(idx)'.*100];

% Correct residual
TGSD(:,2) = TGSD(:,2) + (100-sum(TGSD(:,2)))/size(TGSD,1);

% Write to Fall3D
dlmwrite('fall3d.tgsd', size(TGSD,1), 'delimiter', '\t');

dens = 2500;
spher = 0.9;

F3D = [2.^(TGSD(:,1)), ones(size(TGSD(:,1))).*dens, ones(size(TGSD(:,1))).*spher, TGSD(:,2)./100];
dlmwrite('fall3d.tgsd', F3D, 'delimiter', '\t', '-append');