function plotWind(file, varargin)
% PLOTWIND(path_to_file, ...)
% The function PLOTWIND plots wind profiles saved in 3-columns,
% tab-delimited ascii files containing altitude (m asl), wind speed (m/s) and wind
% direction (degrees from north). 
%
% path_to_file is either a string containing the path to the wind file or a
% cell array containing the individual path to multiple files. If no
% argument is give, a GUI opens.
%
% Name-value pairs of optional input arguments:
% - 'legend':   A cell array of string of the same size as the file paths 
%               containing the legend labels 
%
% Single input flags
% - '-nolegend':  Kills legend


% Define default values
noLegend = 0;
leg = [];

% Check input parameters
if nargin == 0
    [flName,dirName] = uigetfile('*.*', 'Select the wind profile', 'MultiSelect', 'on');
%     if (~ischar(flName) | ~iscell(flName)) & flName == 0; return; end
    if ischar(flName)
        file = {fullfile(dirName, flName)};
    elseif iscell(flName)
        file = fullfile(dirName, flName);
        leg = flName;
    elseif flName == 0
        return;
    end
    
else
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'legend')
            leg = varargin{i+1};
        end
    end
    if nnz(strcmpi(varargin, '-nolegend'))
        noLegend = 1;
    end
end


% Load and store data
if ~iscell(file)
    file = {file};
end
data = cell(length(file), 1);
for i = 1:length(file)
    data{i} = load(file{i});
end

% Prepare figure
figure;
ax1 = subplot(1,2,1);
ax1.Box = 'on';
xlabel(ax1, 'Wind speed (m/s)')
ylabel(ax1, 'Altitude (km)')
hold on

ax2 = subplot(1,2,2);
ax2.Box = 'on';
ax2.XTick = [0, 90, 180, 270];
ax2.XTickLabel = {'N', 'E', 'S', 'W'};
xlabel('Wind direction')
hold on

for iW = 1:length(file)
    plot(ax1, data{iW}(:,2), data{iW}(:,1)./1e3, 'LineWidth', 1);
    plot(ax2, data{iW}(:,3), data{iW}(:,1)./1e3, 'LineWidth', 1);
end

% Check legend
if length(file)>1 & noLegend < 1
    if isempty(leg)
        leg = file;
    end
    legend(ax1, leg, 'location', 'southeast');
end
