% function plotBetaPlume(alpha, beta, height)
% Plot the probability of mass release from the plume following a Beta
% distribution.
% - Alpha and Beta can be entered as single values or as vectors.
% - The height can be entered either as [min max] or as max, in which case
%   the beta distribution is calculated from ground level.
% 
% Examples
%   plotBetaPlume(3, 2, 10)
%   plotBetaPlume([3,3,3],[1,2,3],[1,10])
%
% Written by S. Biass, Jan 2018

function plotBetaPlume(alpha, beta, height, varargin)

if numel(height) == 1
    height = [0, height];
end

if numel(alpha) ~= numel(beta)
    error('Alpha and beta should have the same size')
end

% Control if figure is actually displayed. Useful for inversion
vis = 'on';
if nargin>3 && nnz(strcmpi(varargin, '-noplot'))
    vis = 'off';
end


x = linspace(0,1,25);
h = height(1) + x.*(height(2)-height(1));

stor = zeros(numel(x), numel(alpha), 2);
leg  = cell(numel(alpha),1);
for i = 1:numel(alpha)
    stor(:,i,1) = betapdf(x,alpha(i), beta(i))';
    %stor(:,i,1) = stor(:,i,1)./sum(stor(:,i,1));
    stor(:,i,2) = betacdf(x,alpha(i), beta(i))';
    leg{i}      = ['\alpha = ', num2str(alpha(i)), ', \beta = ', num2str(beta(i))];
end
stor(isinf(stor(:,:,1))) = nan;
stor(:,:,1) = stor(:,:,1)./nansum(stor(:,:,1));

figure('Visible', vis);
subplot(2,1,1)
plot(stor(:,:,1).*100, h);
grid on
xlabel('Release PDF (%)');
ylabel('Plume height')
legend(leg)
grid on

subplot(2,1,2)
plot(stor(:,:,2).*100, h);
grid on
xlabel('Release CDF (%)');
ylabel('Plume height')
grid on