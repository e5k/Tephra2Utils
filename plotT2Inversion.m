function plotT2Inversion(varargin)

if nargin == 0
    [fl,pth] = uigetfile('*.mat');
    data = load(fullfile(pth,fl), 'inversion');
    data = data.inversion;
    ttl = fl;
elseif nargin == 1 && ischar(varargin{1})
    data = load(varargin{1}, 'inversion');
    data = data.inversion;
    [~,fl,ext] = fileparts(varargin{1});
    ttl = strcat(fl,ext);
else
    data = varargin{1};
    ttl = [];
end

h1 = figure(...
    'Units',get(0,'defaultfigureUnits'),...
    'Position',[644 212 900 600],...
    'Visible',get(0,'defaultfigureVisible'),...
    'Color',get(0,'defaultfigureColor'),...
    'CurrentAxesMode','manual',...
    'MenuBar','none',...
    'Name',['Inversion: ', ttl],...
    'NumberTitle','off',...
    'Tag','figure1',...
    'Resize','off');


h2 = uipanel(...
    'Parent',h1,...
    'Title','Plot Tephra2 Inversion',...
    'Units', 'normalized',...
    'Position',[.02 .02 .96 .96] );

h3 = axes(...
    'Parent', h2,...
    'Units', 'normalized',...
    'Position',[.35 .1 .6 .8],...
    'Tag', 'ax',...
    'NextPlot', 'replacechildren');

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .85 .2 .05],...
    'String', 'X Axis:',...
    'HorizontalAlignment','left');

xAx = uicontrol(...
    'Parent', h2,...
    'Style','popupmenu',...
    'Units', 'normalized',...
    'Position',[.02 .8 .2 .05],...
    'Tag', 'xAx',...
    'String', {'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi'},...
    'Value', 2,...
    'Callback', @prepare_data);

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .75 .2 .05],...
    'String', 'Y Axis:',...
    'HorizontalAlignment','left');

yAx = uicontrol(...
    'Parent', h2,...
    'Style','popupmenu',...
    'Units', 'normalized',...
    'Position',[.02 .7 .2 .05],...
    'Tag', 'yAx',...
    'String', {'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi'},...
    'Value', 1,...
    'Callback', @prepare_data);

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .63 .2 .05],...
    'String', 'Minimum fit value:',...
    'HorizontalAlignment','left');

cMin = uicontrol(...
    'Parent', h2,...
    'Style','edit',...
    'Units', 'normalized',...
    'Position',[.02 .6 .2 .05],...
    'String', 'aa',...
    'Tag', 'cMin',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .53 .2 .05],...
    'String', 'Minimum fit value:',...
    'HorizontalAlignment','left');

cMax = uicontrol(...
    'Parent', h2,...
    'Style','edit',...
    'Units', 'normalized',...
    'Position',[.02 .5 .2 .05],...
    'String', 'aa',...
    'Tag', 'cMax',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .43 .2 .05],...
    'String', 'X Steps:',...
    'HorizontalAlignment','left');

xStep = uicontrol(...
    'Parent', h2,...
    'Style','edit',...
    'Units', 'normalized',...
    'Position',[.02 .4 .2 .05],...
    'String', '20',...
    'Tag', 'xStep',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

uicontrol(...
    'Parent', h2,...
    'Style','text',...
    'Units', 'normalized',...
    'Position',[.02 .33 .2 .05],...
    'String', 'Y Steps:',...
    'HorizontalAlignment','left');

yStep = uicontrol(...
    'Parent', h2,...
    'Style','edit',...
    'Units', 'normalized',...
    'Position',[.02 .3 .2 .05],...
    'String', '20',...
    'Tag', 'yStep',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

update = uicontrol(...
    'Parent', h2,...
    'Style','pushbutton',...
    'Units', 'normalized',...
    'Position',[.02 .05 .1 .05],...
    'String', 'Update',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

export = uicontrol(...
    'Parent', h2,...
    'Style','pushbutton',...
    'Units', 'normalized',...
    'Position',[.13 .05 .1 .05],...
    'String', 'Export',...
    'HorizontalAlignment','left',...
    'Callback', @prepare_data);

h1.UserData = data;
cMin.String = num2str(min(data.fitTable.Fit));
cMax.String = num2str(max(data.fitTable.Fit));

prepare_data(update)

function prepare_data(hObject, ~, ~)

fig     = hObject.Parent.Parent;
data    = fig.UserData;

% Retrieve data
xAx     = findobj(fig, 'Tag', 'xAx');
xAx     = xAx.String{xAx.Value};

yAx     = findobj(fig, 'Tag', 'yAx');
yAx     = yAx.String{yAx.Value};

xStep   = findobj(fig, 'Tag', 'xStep');
xStep   = str2double(xStep.String);

yStep   = findobj(fig, 'Tag', 'yStep');
yStep   = str2double(yStep.String);

cMin    = findobj(fig, 'Tag', 'cMin');
cMin    = str2double(cMin.String);

cMax    = findobj(fig, 'Tag', 'cMax');
cMax    = str2double(cMax .String);


% Scale parameters
if strcmp(xAx, 'Mass') || strcmp(yAx, 'Mass')
    data.fitTable.Mass = log10(data.fitTable.Mass);
end

if strcmp(xAx, 'Height') || strcmp(yAx, 'Height')
    data.fitTable.Height = data.fitTable.Height./1e3;
end

% Interpolate
x = linspace(min(data.fitTable.(xAx)), max(data.fitTable.(xAx)), xStep);
% x = [x(1)-(x(2)-x(1)), x, x(end)+(x(end)-x(end-1))];
y = linspace(min(data.fitTable.(yAx)), max(data.fitTable.(yAx)), yStep);
% y = [y(1)-(y(2)-y(1)), y, y(end)+(y(end)-y(end-1))];

[xi,yi] = meshgrid(x,y) ;
% zi      = griddata(data.fitTable.(xAx),data.fitTable.(yAx),data.fitTable.Fit,xi,yi, 'cubic') ;

F = scatteredInterpolant(data.fitTable.(xAx),data.fitTable.(yAx),data.fitTable.Fit, 'natural');
zi = F(xi,yi);

if strcmp(hObject.String, 'Update') | strcmp(hObject.Tag, 'xAx') | strcmp(hObject.Tag, 'yAx')...
         | strcmp(hObject.Tag, 'cMin') | strcmp(hObject.Tag, 'cMax') | strcmp(hObject.Tag, 'xStep') | strcmp(hObject.Tag, 'yStep')
    ax = findobj(fig, 'Type', 'Axes');
else
    f = figure;
    ax = axes('Parent', f);
end

cla(ax, 'reset')

imagesc(ax, x, y, zi);
hold on
plot(data.fitTable.(xAx),data.fitTable.(yAx), 'xk')

if numel(data.fitTable.(yAx))>=10
    bnd = 10;
else
    bnd = numel(data.fitTable.(yAx));
end

plot(data.fitTable.(xAx)(1:bnd),data.fitTable.(yAx)(1:bnd), '.r')
text(data.fitTable.(xAx)(1:bnd),data.fitTable.(yAx)(1:bnd), cellstr(sprintfc('%d',(1:bnd)')), 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');

axis square
ax.YDir = 'Normal';
xlabel(xAx);
ylabel(yAx);

c = colorbar;
set(get(c,'ylabel'),'string','Fit','fontsize',10);
caxis([cMin, cMax])

    
