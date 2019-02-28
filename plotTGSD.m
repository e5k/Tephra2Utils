function plotTGSD(mu, sigma)

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
