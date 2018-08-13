function data = runT2(varargin)
    
run.name    = 'test';
run.mode    = 0;                % 0: Deterministic
                                % 1: Probabilistic
run.mass    = 0;                % 0: Mass if specified
                                % 1: Mass computed from height, wind and duration
                                % 2: Mass computed from thermal theory - to add

grid.pth    = 'gr';%{'1.grd','2.grd'};               %                                 
grid.xLim   = [-50000 250000];  % X limits (UTM, m)
grid.yLim   = [-75000 75000];   % Y limits (UTM, m)
grid.res    = [1000,2000];           % Grid resolution (m)
grid.elev   = [0,100];              % Mean elevation
grid.zone   = [];               % Zone

vent.x      = 0;
vent.y      = 0;
vent.z      = 1;

wind.pth    = [];               % String, cell array of strings, path
wind.dir    = 90;               % Degree from North, wind direction
wind.trop   = 12000;            % Tropopause height (m asl)
wind.strat  = 15000;            % Stratosphere height (m asl)
wind.vel    = 5:5:15;           % Velocity at tropopause (m/s)
wind.model  = 1;                % 1: WM1 of Bonadonna and Phillips
                                % 2: WM2 of Bonadonna and Phillips

tgsd.pth    = [];
tgsd.med    = -1;
tgsd.std    = 2;
tgsd.agg    = 0.5;
tgsd.minD   = 5;

conf.pth           = [];
conf.plumeHeight   = 5000:5000:20000;
conf.mass          = logspace(11,13,10);% 1e11:.5e11:1e13;
conf.duration      = 1; % Hours
conf.eddyConst     = 0.04;
conf.diffCoef      = 112;
conf.ftt           = 1563;
conf.lithDens      = 2600;
conf.pumDens       = 513;
conf.colStep       = 50;
conf.partStep      = 50;
conf.plumeAlpha    = 3;
conf.plumeBeta     = 2;

% Check input structure
if nargin > 0 && isstruct(varargin{1})
    data    = varargin{1};
    run     = data.run;
    grid    = data.grid;
    vent    = data.vent;
    wind    = data.wind;
    tgsd    = data.tgsd;
    conf    = data.conf;
end
    
% Check compilation
if nnz(strcmp(varargin,'-compile')) > 0
    compile = 1;
else
    compile = 0;
end

%% Prepare output
mkdir(run.name);
mkdir(fullfile(run.name, 'CONF'));
mkdir(fullfile(run.name, 'GRID'));
mkdir(fullfile(run.name, 'WIND'));
mkdir(fullfile(run.name, 'OUT'));
mkdir(fullfile(run.name, 'TGSD'));

elements    = cell(30,1);
elementsN   = cell(30,1);
elemC       = 1;

%% Grid
elemS   = elemC; % Start counter
% In case the path to the grid files is specified
if ~isempty(grid.pth)
    grid.elem    = {'pth'};
    grid.elemF   = {'Grid name'};
    
    % Check if the path is a directory
    if isdir(grid.pth)
        grid.pth = dir2pth(grid.pth);  
    end
    
    [elements, elementsN, elemC, grid] ...
        = prepareInput(elements, elementsN, elemS, elemC, grid, 'grid', '.grd',1);
% Else get combinations and create grids
else
    grid.elem    = {'res', 'elev'};                          % Variable name
    grid.elemF   = {'Grid resolution', 'Grid elevation'};    % Full name
    [elements, elementsN, elemC, grid] ...
        = prepareInput(elements, elementsN, elemS, elemC, grid, 'grid', '.grd',0);
    
    % Create grids
    for i = 1:length(grid.P.res)
    makeGrid(grid.xLim(1):grid.P.res(i):grid.xLim(2),...
        grid.yLim(1):grid.P.res(i):grid.yLim(2),...
        grid.P.elev(i),...
        fullfile(run.name, 'GRID', grid.N{i}));
    end
end
                                         




%% Wind
if isempty(wind.pth)
    wind.elem    = {'dir', 'vel', 'trop', 'strat', 'model'}; % Variable name
    wind.elemF   = {'Wind direction', 'Wind velocity', 'Tropopause height', 'Stratosphere height', 'Wind model'};
    elemS   = elemC; % Start counter

    [elements, elementsN, elemC, wind] ...
            = prepareInput(elements, elementsN, elemS, elemC, wind, 'wind', '.wnd',0);
    
    for i = 1:length(wind.P.dir)
        makeWindProfile(wind.P.dir(i), wind.P.vel(i), wind.P.trop(i), wind.P.strat(i), wind.P.model(i), fullfile(run.name, 'WIND', wind.N{i}));
    end
else
    
end
    
%% TGSD
if isempty(tgsd.pth)
    tgsd.elem    = {'med', 'std', 'agg', 'minD'};
    tgsd.elemF   = {'Median TGSD', 'Std TGSD', 'Aggregation coef.', 'Min aggregate diam'};
    elemS   = elemC;

    [elements, elementsN, elemC, tgsd] ...
            = prepareInput(elements, elementsN, elemS, elemC, tgsd, 'tgsd', '.gsd',0);
    
    for i = 1:length(tgsd.P.med)
        makeTGSD(tgsd.P.med(i), tgsd.P.std(i), tgsd.P.agg(i), tgsd.P.minD(i), fullfile(run.name, 'TGSD', tgsd.N{i}));
    end
end

%% Configuration file
conf.windVel    = wind.vel; 
conf.windTrop   = wind.trop; 

var     = {'plumeHeight', 'mass', 'duration', 'eddyConst', 'diffCoef', 'ftt', 'lithDens', 'pumDens', 'colStep', 'partStep', 'plumeAlpha', 'plumeBeta', 'windVel','windTrop'};
varF    = {'Plume height', 'Mass', 'Duration', 'Eddy constant', 'Diffusion coef', 'FTT', 'Lithic density', 'Pumice density', 'Column steps', 'Particle steps', 'Alpha', 'Beta', 'wind', 'trop'};

% Index used to define what variables are associated for a given run.mass
if run.mass == 0
    idxE         = ~strcmp(var, 'duration');
elseif run.mass == 1
    idxE         = ~strcmp(var, 'mass');
end

conf.elem    = var(idxE);
conf.elemF   = varF(idxE);
elemS        = elemC;

[~,~,~, conf] ...
    = prepareInput(elements, elementsN, elemS, elemC, conf, 'conf', '.conf',0);

% Calculate MER
MER = zeros(length(conf.P.plumeHeight),1);
for i = 1:length(conf.P.plumeHeight)
    MER(i) = get_MER_DB12(conf.P.plumeHeight(i)-vent.z, ...
        conf.P.windVel(i),...
        conf.P.windTrop(i));
end
% Add MER to the conf structure
conf.P.MER  = MER;
conf.elem   = [conf.elem, 'MER'];
conf.elemF  = [conf.elemF, 'MER'];

% Index used to define what variables are associated for a given run.mass
if run.mass == 0
    conf.P.duration = conf.P.mass./conf.P.MER./3600;
    conf.elem   = [conf.elem, 'duration'];
    conf.elemF  = [conf.elemF, 'Duration'];
elseif run.mass == 1
    conf.P.mass = conf.P.MER .* conf.P.duration .* 3600;
    conf.elem   = [conf.elem, 'mass'];
    conf.elemF  = [conf.elemF, 'Mass'];
end

makeConf(conf, vent, run.name);

%% Run the hard way
sConf = size(conf.N, 1);
sGrid = size(grid.N, 1);
sWind = size(wind.N, 1);
sTgsd = size(tgsd.N, 1);

% This vector contains the indices to conf, grid, wind and tgsd in each
% column
idxStor = zeros(sConf * sGrid * sWind * sTgsd,4);
T2Stor  = cell(sConf * sGrid * sWind * sTgsd,1);
count   = 1;

% Create a master structure containing all elements
masterRun = struct;
for i = 1:length(conf.elem); masterRun.(conf.elem{i}) = []; end
for i = 1:length(grid.elem); masterRun.(grid.elem{i}) = []; end
for i = 1:length(wind.elem); masterRun.(wind.elem{i}) = []; end
for i = 1:length(tgsd.elem); masterRun.(tgsd.elem{i}) = []; end

for iConf = 1:sConf
    for iGrid = 1:sGrid
        for iWind = 1:sWind
            for iTgsd = 1:sTgsd
                
                for i = 1:length(conf.elem); masterRun.(conf.elem{i})(count) = conf.P.(conf.elem{i})(iConf); end
                for i = 1:length(grid.elem); masterRun.(grid.elem{i})(count) = grid.P.(grid.elem{i})(iGrid); end
                for i = 1:length(wind.elem); masterRun.(wind.elem{i})(count) = wind.P.(wind.elem{i})(iWind); end
                for i = 1:length(tgsd.elem); masterRun.(tgsd.elem{i})(count) = tgsd.P.(tgsd.elem{i})(iTgsd); end
                
                idxStor(count,:) = [iConf, iGrid iWind, iTgsd];
                
                % Adapts file path
                if isempty(grid.pth); gridP = [run.name, filesep, 'GRID', filesep]; else; gridP = []; end
                if isempty(wind.pth); windP = [run.name, filesep, 'WIND', filesep]; else; windP = []; end
                if isempty(tgsd.pth); tgsdP = [run.name, filesep, 'TGSD', filesep]; else; tgsdP = []; end
                if isempty(conf.pth); confP = [run.name, filesep, 'CONF', filesep]; else; confP = []; end

                T2Stor{count}    = ['tephra2-2012 ', confP, conf.N{iConf}, ' ', gridP, grid.N{iGrid}, ' ', windP, wind.N{iWind}, ' ', tgsdP, tgsd.N{iTgsd}, ' > ', run.name, filesep, 'OUT', filesep, num2str(count, '%5.0f'), '.out'];
                
                % Appends ./ on unix
                if ~ispc
                    T2Stor{count}    = ['./', T2Stor{count}];
                end
                count = count+1;
            end
        end
    end
end

if compile == 1
    compileT2;
end

data.run    = run;
data.grid   = grid;
data.vent   = vent;
data.wind   = wind;
data.tgsd   = tgsd;
data.conf   = conf;
data.master = masterRun;
data = gatherParam(data,T2Stor);

save([run.name, '.mat'], 'data', 'idxStor')

runIt(T2Stor);


% Cleanup
elements  = elements(~cellfun(@isempty, elements));
elementsN = elementsN(~cellfun(@isempty, elements));


function data = gatherParam(data,T2Stor)

param   = [data.conf.elem, data.wind.elem, data.grid.elem, data.tgsd.elem];
%stor    = zeros(size(T2Stor,1), numel(param));

for i = 1:size(T2Stor,1)
    st = strsplit(T2Stor{i}, '_');
    for j = 1:length(param)
        for k = 1:length(st)
            if ~isempty(regexp(st{k}, param{j}, 'once'))
                tmp = str2double(regexp(st{k},'\d+\.?\d*', 'match'));
                stor.(param{j})(i) = tmp(1);
               %stor(i,j) = str2double(regexp(st{k},'\d+\.?\d*', 'match'));
            end
        end
    end
end
data.runlist = stor;


% For each of grid, conf, gs and wind, return the all the combinations of
% variables and generate file names
% type defines if filenames are input (1) or values (0)
function  [elements, elementsN, elemC, struc] = prepareInput(elements, elementsN, elemS, elemC, struc, fl, ext, type)
for i = 1:length(struc.elem)
    elements{elemC}     = struc.(struc.elem{i});
    elementsN{elemC}    = struc.elemF{i};
    elemC               = elemC+1;
end

% Create combinations
if type == 0
    struc.Pt  = makePossibilities(elements(elemS:elemC-1));   
    struc.N   = cellfun(@num2str, num2cell(struc.Pt), 'UniformOutput', false);
    struc.N   = makeName(struc.N, struc.elem, fl, ext);
% Else use file path
else
    struc.Pt = elements{elemS:elemC-1};
    struc.N  = elements{elemS:elemC-1};
    struc.Pt = (1:numel(struc.Pt))';
end
% Here convert all variables into structure fields
for i = 1:length(struc.elem)
    struc.P.(struc.elem{i}) = struc.Pt(:,i);
end

% Create file names
function outN = makeName(N, elem, fl, ext)
for i = 1:size(N,2)
    N(:,i) = strcat(['_', elem{i}], '-', N(:,i));
end
outN = cell(size(N,1),1);
for i = 1:size(N,1)
    outN{i} = strcat(fl, N{i,:}, ext);
end

function result = makePossibilities(elements)
combinations        = cell(1, numel(elements)); %set up the varargout result
[combinations{:}]   = ndgrid(elements{:});
combinations        = cellfun(@(x) x(:), combinations,'uniformoutput',false); 
result              = [combinations{:}]; % NumberOfCombinations by N matrix. Each row is uni

% % Generation of input parameters
% Create grid
function makeGrid(xVec,yVec,elev,name)
[X,Y] = meshgrid(xVec, yVec);
utm   = [reshape(X, numel(X),1), reshape(Y, numel(X),1), reshape(ones(numel(X),1).*elev, numel(X),1)];
dlmwrite(name, utm, 'delimiter', '\t')

% Create wind profiles
function makeWindProfile(dir, speed, trop, strat, model, name)
dummyH      = [0:50:500, 600:100:2000, 2500:500:10000, 11000:1000:25000, 27500:2500:40000];
dummyH(1)   = 1; % Remove a 0 elevation
vel         = zeros(size(dummyH));

[~, idxT]   = min(abs(dummyH-trop));
[~, idxS]   = min(abs(dummyH-strat));
[~, idx20]  = min(abs(dummyH-20000));

vel(1:idxT) = interp1([0,trop], [0,speed], dummyH(1:idxT));

if model == 1
    vel(idxT+1:idxS) = interp1([trop,strat], [speed, .75*speed], dummyH(idxT+1:idxS));
    vel(idxS+1:end)  = ones(size(vel(idxS+1:end))) .* .75*speed;
else    
    vel(idxT+1:idx20) = interp1([trop,20000], [speed, .1*speed], dummyH(idxT+1:idx20));
    vel(idx20+1:end)  = ones(size(vel(idx20+1:end))) .* .1*speed;
end
dlmwrite(name, [dummyH', vel', ones(size(vel')).*dir], 'delimiter', '\t');

% Create TGSD
function makeTGSD(med, std, agg, maxD, name)

diam    = -15:15;
i1      = diam >= -1 & diam < maxD;
i2      = diam >= maxD & diam <= diam(end);

wt      = normpdf(diam, med, std);
%wt      = round(wt,4);

wt(i2)  = wt(i2) - wt(i2)*agg;
wt(i1)  = wt(i1) + sum(wt(i2).*agg)/nnz(i1);
wt      = round(wt,4);
diam    = diam(wt>0);
wt      = wt(wt>0); 

% Convert to cdf
for i = 2:length(wt)
   wt(i) = wt(i-1)+wt(i); 
end

% Correct rounding
wt(wt==wt(end)) = 1;

dlmwrite(name, [diam', wt'], 'delimiter', '\t', 'precision', 4);

% Write configuration files
function makeConf(conf, vent, name)
for i = 1:size(conf.N,1)
    fid = fopen(fullfile(name, 'CONF', conf.N{i}), 'wt');
    fprintf(fid,...
        'PLUME_HEIGHT\t%d\nERUPTION_MASS\t%d\nVENT_EASTING\t%d\nVENT_NORTHING\t%d\nVENT_ELEVATION\t%d\nEDDY_CONST\t%d\nDIFFUSION_COEFFICIENT\t%d\nFALL_TIME_THRESHOLD\t%d\nLITHIC_DENSITY\t%d\nPUMICE_DENSITY\t%d\nCOL_STEPS\t%d\nPART_STEPS\t%d\nPLUME_MODEL\t%d\nALPHA\t%d\nBETA\t%d\n',...
        conf.P.plumeHeight(i),...
        conf.P.mass(i),...
        vent.x, vent.y, vent.z,...
        conf.P.eddyConst(i),...
        conf.P.diffCoef(i),...
        conf.P.ftt(i),...
        conf.P.lithDens(i),...
        conf.P.pumDens(i),...
        conf.P.colStep(i),...
        conf.P.partStep(i),...
        2,...
        conf.P.plumeAlpha(i),...
        conf.P.plumeBeta(i));
    fclose(fid);
end

% Compile Tephra2
function compileT2
% Retrieve path
pth     = pwd;
mod_pth = [pwd, filesep, 'Tephra2', filesep, 'forward_src', filesep];

% On PCs, add the Cygwin folder to the PATH environment
if ispc
    path1 = getenv('PATH');
    if strcmp(computer('arch'), 'win64')
        pathC = 'C:\cygwin64\bin';
    elseif strcmp(computer('arch'), 'win32')
        pathC = 'C:\cygwin\bin';
    end
        
    if isempty(regexp(path1,pathC,'once'))
        if isdir(pathC)
            setenv('PATH', [path1,';',pathC]);
        else
            choice = questdlg('Did you already install CYGWIN?', ...
                'CYGWIN', ...
                'Yes','No','Yes');
            % Handle response
            switch choice
                case 'Yes'
                    fname = uigetdir('C:\', 'Select the cygwin\bin directory');
                    setenv('PATH', [path1,';',fname]);
                case 'No'
                    url('https://cygwin.com/install.html');
                    return
            end
        end
    end
end
        

% Compiles the model and runs it
cd(mod_pth);                            % Navigates to the makefile
system('make clean');
[stat, cmd_out] = system('make');       % Compiles TEPHRA2

if stat == 0                            % If compilation ok
    cd(pth);
    tmp = dir('Tephra2/tephra2-2012*');
    movefile(['Tephra2/', tmp(1).name], pth);
else
    cd(pth);
    error('There was a problem compiling TEPHRA2...\n%s\n', cmd_out);
end

% Run Tephra2
function runIt(stor)

if license('test','Distrib_Computing_Toolbox')
    %n=str2num(getenv('SLURM_CPUS_PER_TASK'));
    % create a parpool
    %parpool(n);
    %parpool(4)
    parfor i = 1:size(stor,1)
        display(['Run ', num2str(i), '/',num2str(size(stor,1))]);
        system(stor{i});
        home;
    end
else 
    for i = 1:size(stor,1)
        display(['Run ', num2str(i), '/',num2str(size(stor,1))]);
        system(stor{i});
        home;
    end
end
home;
delete('node_');
delete('plume2.dat');
disp('Modelling finished!');

% Converts a directory to a cell array to files
function names = dir2pth(pth)
fl    = dir(pth);
names = {fl.name};
names = names{~[fl.isdir]};