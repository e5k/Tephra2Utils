function [E,N,M,varargout] = plotT2(data, varargin)
% PLOTT2(path_to_file, ...)
% The function PLOTT2 plots the output of Tephra2 on a map
% path_to_file is a string containing the path to the output file
%
% Name-value pairs of optional input arguments:
% - 'zone':     UTM zone of the study area (double). If specified and the 
%               function utm2ll is available, the dataset is converted to 
%               geographic coordinates. If the function plot_google_map is
%               available, the Google Map background is plotted
% - 'plot':     Defines if accumulation is plotted as log10 ('log10,
%               default) or linear ('linear')
% - 'contours': Values of the contours (vector, kg/m2), by default 
%               [.1, 1, 10, 100, 1000]
%
% Outputs
% - E, N, M     Gridded versions of easting, northing and accumulation,
%               respectively
% - If 'zone' is specified, also returns the gridded latitude/longitude
%
% Example
% [E,N,M]           = plotT2('example.out');
% [E,N,M,LAT,LON]   = plotT2('example.out', 'zone', 18);
% [E,N,M]           = plotT2('example.out', 'plot', 'linear');
% [E,N,M]           = plotT2('example.out', 'contours', [10, 50, 100]);

%% Parse optional input arguments
zone    = [];
plt     = 'log';
cnt     = [.1, 1, 10, 100, 1000];
GISname = [];

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'zone')
        zone = varargin{i+1};
        i = i+1;
    elseif strcmpi(varargin{i}, 'plot')
        plt = varargin{i+1};
        i = 1+1;
    elseif strcmpi(varargin{i}, 'contours')
        cnt = varargin{i+1};
        i = 1+1;
    elseif strcmpi(varargin{i}, 'outName')
        GISname = varargin{i+1};
        i = 1+1;
    end
end

%% Calculations

% Load output file
try
    data = load(data);                  % If the Tephra2 output file does not have a header
catch ME
    data = dlmread(data, ' ', 1, 0);
end

% Prepare the data
data    = data(:,1:4);                  % Remove GS data
data    = sortrows(data, [1,2]);        % Sort along northing

% Define grid size
xn      = length(unique(data(:,1)));    % Number of easting values
yn      = length(unique(data(:,2)));    % Number of northing values

% Reshape data into a grid
E       = reshape(data(:,1), yn, xn);
N       = reshape(data(:,2), yn, xn);
M       = reshape(data(:,4), yn, xn);

M(M<cnt(1)) = nan;                      % Remove low values

% If a UTM zone is passed, convert to lat/lon
if ~isempty(zone)
    if exist('utm2ll.m', 'file')
        [LT,LN] = utm2ll(data(:,1), data(:,2), zone);
        LT = reshape(LT, xn, yn);
        LN = reshape(LN, xn, yn);
    else
        warning('The function utm2ll is not available. Download it from https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll');
    end
end

%% Plot
% Define if projected or geographic
if exist('LT', 'var')
    XX  = LN;
    YY  = LT;
    xl  = 'Longitude';
    yl  = 'Latitude';
else 
    XX  = E;
    YY  = N;
    xl  = 'Easting';
    yl  = 'Northing';   
end

% Define if data is plotted as log
if strcmpi(plt, 'log')
    M   = log10(M);
    cnt = log10(cnt);
    zl  = 'Log_{10} accumulation (kg/m^2)';
else
    zl  = 'Tephra accumulation (kg/m^2)';
end

% Plot figure
res         = (XX(1,2)-XX(1,1))/2;
figure;
hd          = pcolor(XX-res,YY-res,M); shading flat; hold on;
[c,h]       = contour(XX,YY,M,cnt, 'Color', 'k');
clabel(c,h, cnt, 'LabelSpacing', 1000, 'FontWeight', 'bold')
set(hd, 'FaceAlpha', 0.5)

% 
if exist('plot_google_map', 'file') && exist('LT', 'var')
    plot_google_map('maptype', 'terrain');
    
    % Plot grid extent
    gX = [XX(1,1), XX(1,end), XX(end,end), XX(end,1), XX(1,1)];
    gY = [YY(1,1), YY(1,end), YY(end,end), YY(end,1), YY(1,1)];
    plot(gX, gY, '-r', 'linewidth',0.5);
    
elseif ~exist('plot_google_map', 'file') && exist('LT', 'var')
    warning('The function plot_google_map is not available. Download it from https://uk.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot-google-map');
else
    axis equal tight;
end

c = colorbar;

% Labels
ylabel(c, zl);
xlabel(xl);
ylabel(yl);

set(gca, 'Layer', 'top');

%% Prepare output
% If lat/lon, output as arguments 4 and 5
if exist('LT', 'var')
    varargout{1} = LT;
    varargout{2} = LN;
end
