function [E,N,M,Carea,Pmass,cont,varargout] = plotT2(data, varargin)
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
% - 'vent':     Vent coordinates, entered as [easting, northing]
% - 'minVal':   Minimum value to be represented on the continuous color surface (kg/m2)
% - 'points':   Specify points to plot as a 2-columns [easting,northing] matrix. If entered as a 3-columns [easting,northing,value] matrix, value is used as a lable
%
% Single input flags
% - 'noplot':   Does not display a map
%
% Outputs
% [E,N,M,A]         = plotT2(...)
%               Return gridded easting, northing, mass accumulation and isomass areas (km2)
% [E,N,M,A,Lat,Lon] = plotT2(...)
%               If 'zone' is specified, also returns the gridded latitude/longitude
%
% Examples
% [E,N,M,A]         = plotT2('example.out');
% [E,N,M,A,LAT,LON] = plotT2('example.out', 'zone', 18);
% [E,N,M,A]         = plotT2('example.out', 'plot', 'linear');
% [E,N,M,A]         = plotT2('example.out', 'contours', [10, 50, 100]);

%% Parse optional input arguments
zone    = [];
plt     = 'log';
cnt     = [.1, 1, 10, 100, 1000];
noPlot  = 0;
vent    = [];
points  = [];
minVal  = 0;

if nargin == 0
    [flName,dirName] = uigetfile('*.*', 'Select the Tephra2 output file');
    if flName == 0; return; end
    data = fullfile(dirName, flName);
else
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
        elseif strcmpi(varargin{i}, 'vent')
            vent = varargin{i+1};
            i = 1+1;
        elseif strcmpi(varargin{i}, 'minVal')
            minVal = varargin{i+1};
            i = 1+1;
        elseif strcmpi(varargin{i}, 'points')
            points = varargin{i+1};
            i = 1+1;
        end
    end
    if nnz(strcmpi(varargin, '-noplot'))
        noPlot = 1;
    end
end

% Temporarly remove warnings
warning off backtrace % turn all warnings off

%% Calculations
% Load output file
flName = data;
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

% If a UTM zone is passed, convert to lat/lon
if ~isempty(zone)
    if exist('utm2ll.m', 'file')
        [LT,LN]     = utm2ll(data(:,1), data(:,2), zone); 
        if ~isempty(vent); [LTv,LNv]     = utm2ll(vent(:,1), vent(:,2), zone); end
        if ~isempty(points); [LTp,LNp]   = utm2ll(points(:,1), points(:,2), zone); end
        LT = reshape(LT, yn, xn);
        LN = reshape(LN, yn, xn);
    else
        warning('The function utm2ll is not available. Download it from https://www.mathworks.com/matlabcentral/fileexchange/45699-ll2utm-and-utm2ll');
    end
end

% Retrieve the area and the mass at each isomass
[Carea, Pmass, cont] = getArea(E(1,:), M, cnt);


%% Plot
% Define if projected or geographic
if exist('LT', 'var')
    XX  = LN;
    YY  = LT;
    xl  = 'Longitude';
    yl  = 'Latitude';
    if ~isempty(vent);   V   = [LNv, LTv]; end
    if ~isempty(points); P   = [LNp, LTp]; end
else 
    XX  = E;
    YY  = N;
    xl  = 'Easting';
    yl  = 'Northing';   
    if ~isempty(vent);   V   = vent; end    
    if ~isempty(points); P   = points(:,1:2); end
end

% Remove values smaller than minVal
M(M<minVal) = nan;

% Define if data is plotted as log
if strcmpi(plt, 'log')
    M   = log10(M);
    cnt = log10(cnt);
    zl  = 'Log_{10} accumulation (kg/m^2)';
else
    zl  = 'Tephra accumulation (kg/m^2)';
end

% Plot figure
if noPlot == 0
    figure;
    res         = (XX(1,2)-XX(1,1))/2;
    hd          = pcolor(XX-res,YY-res,M); shading flat; hold on;
    [c,h]       = contour(XX,YY,M,cnt, 'Color', 'k');
    clabel(c,h, cnt, 'LabelSpacing', 1000, 'FontWeight', 'bold')
    set(hd, 'FaceAlpha', 0.5)
    % Plot vent
    if ~isempty(vent)
        plot(V(1), V(2), '^k', 'MarkerFaceColor', 'r', 'MarkerSize', 10)
    end
    % Plot points
    if ~isempty(points)
        plot(P(:,1), P(:,2), 'ok', 'MarkerFaceColor', 'w', 'MarkerSize', 5);
        if size(points,2) == 3
            text(P(:,1), P(:,2), cellstr(num2str(points(:,3), '%.1f')));
        end
    end
    
    if exist('plot_google_map', 'file') && exist('LT', 'var')
        plot_google_map('maptype', 'terrain');
        % Plot grid extent
        gX = [XX(1,1), XX(1,end), XX(end,end), XX(end,1), XX(1,1)];
        gY = [YY(1,1), YY(1,end), YY(end,end), YY(end,1), YY(1,1)];
        plot(gX-res, gY-res, '-r', 'linewidth',0.5);
    elseif ~exist('plot_google_map', 'file') && exist('LT', 'var')
        warning('The function plot_google_map is not available. Download it from https://uk.mathworks.com/matlabcentral/fileexchange/27627-zoharby-plot-google-map');
    else
        axis equal tight;
    end

    c = colorbar;

    % Labels
    title(flName,'Interpreter', 'none');
    ylabel(c, zl);
    xlabel(xl);
    ylabel(yl);
    set(gca, 'Layer', 'top');
end

%% Prepare output
if nargout > 0
    varargout = cell(nargout,1);
    for i = 1:nargout
        if      i == 1; varargout{i} = E;
        elseif  i == 2; varargout{i} = N;
        elseif  i == 3; varargout{i} = M;
        elseif  i == 4 && exist('LT', 'var'); varargout{i} = LT;
        elseif  i == 5 && exist('LT', 'var'); varargout{i} = LN;
        end
    end
end
warning on backtrace % turn all warnings on

% Calculates the area for all contours
function [Carea, Pmass, cont] = getArea(E, M, cnt)
 
res = E(1,2)-E(1,1);
Carea = zeros(numel(cnt),1);
Pmass = zeros(numel(cnt),1);

for i = 1:numel(cnt)
    idx = M>=cnt(i);
    Carea(i) = (nnz(idx).*res^2)/1e6;
    Pmass(i) = sum(M(idx)).*res^2;
    % Check that a given contour isn't truncated
    if nnz(idx(1,:))>0 | nnz(idx(end,:))>0 | nnz(idx(:,1))>0 | nnz(idx(:,end))>0 
        Carea(i) = nan;
        Pmass(i) = nan;
        warning('The contour %d extends outside of the domain', cnt(i));
    end
    % Check if the contour exists
    if nnz(idx) == 0
        Carea(i) = nan;
        Pmass(i) = nan;
        warning('The contour %d does not exist', cnt(i));
%     else
%         Carea(i) = (nnz(idx).*res^2)/1e6;
%         Pmass(i) = sum(M(idx)).*res^2;
    end
end

idx = isnan(Carea);
Carea(idx) = [];
Pmass(idx) = [];
cont       = cnt(~idx);