function runT2(varargin)

run.name    = 'test';
run.mode    = 0;                % 0: Deterministic
                                % 1: Probabilistic
run.mass    = 0;                % 0: Mass if specified
                                % 1: Mass computed from height, wind and duration
                                % 2: Mass computed from thermal theory - to add
                                
grid.xLim   = [-50000 250000];  % X limits (UTM, m)
grid.yLim   = [-75000 75000];   % Y limits (UTM, m)
grid.res    = [1000,2000];             % Grid resolution (m)
grid.elev   = [0,10];                % Mean elevation
grid.zone   = [];               % Zone

vent.x      = 0;
vent.y      = 0;
vent.z      = 1;

wind.pth    = [];               % If not empty, then use it
wind.dir    = 90;               % Degree from North, wind direction
wind.trop   = 12000;            % Tropopause height (m asl)
wind.strat  = 15000;            % Stratosphere height (m asl)
wind.vel    = [5:5:30];         % Velocity at tropopause (m/s)
wind.model  = 1;                % 1: WM1 of Bonadonna and Phillips
                                % 2: WM2 of Bonadonna and Phillips

tgsd.pth    = [];
tgsd.med    = -1;
tgsd.std    = 2;
tgsd.agg    = 0.5;
tgsd.minD   = 5;


conf.plumeHeight   = [1000, 2500, 5000:5000:30000];
conf.mass          = 1e10;
conf.duration      = 1; % Hours
conf.eddyConst     = 0.04;
conf.diffCoef      = 112;
conf.ftt           = 1563;
conf.lithDens      = [2600,3000];
conf.pumDens       = 513;
conf.colStep       = 50;
conf.partStep      = 50;
conf.plumeAlpha    = 3;
conf.plumeBeta     = 2;


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
grid.elem    = {'res', 'elev'};                          % Variable name
grid.elemF   = {'Grid resolution', 'Grid elevation'};    % Full name
elemS   = elemC;                                    % Start counter

[elements, elementsN, elemC, grid] ...
        = prepareInput(elements, elementsN, elemS, elemC, grid, 'grid', '.grd');

for i = 1:size(grid.P,1)
    makeGrid(grid.xLim(1):grid.P(i,1):grid.xLim(2),...
        grid.yLim(1):grid.P(i,1):grid.yLim(2),...
        grid.P(i,2),...
        fullfile(run.name, 'GRID', grid.N{i}));
end

%% Wind
if isempty(wind.pth)
    wind.elem    = {'dir', 'vel', 'trop', 'strat', 'model'}; % Variable name
    wind.elemF   = {'Wind direction', 'Wind velocity', 'Tropopause height', 'Stratosphere height', 'Wind model'};
    elemS   = elemC; % Start counter

    [elements, elementsN, elemC, wind] ...
            = prepareInput(elements, elementsN, elemS, elemC, wind, 'wind', '.wnd');
    
    for i = 1:size(wind.P,1)
        makeWindProfile(wind.P(i,1), wind.P(i,2), wind.P(i,3), wind.P(i,4), wind.P(i,5), fullfile(run.name, 'WIND', wind.N{i}));
    end
end
    
%% TGSD
if isempty(tgsd.pth)
    tgsd.elem    = {'med', 'std', 'agg', 'minD'};
    tgsd.elemF   = {'Median TGSD', 'Std TGSD', 'Aggregation coef.', 'Min aggregate diam'};
    elemS   = elemC;

    [elements, elementsN, elemC, tgsd] ...
            = prepareInput(elements, elementsN, elemS, elemC, tgsd, 'tgsd', '.gsd');
    
    for i = 1:size(tgsd.P,1)
        makeTGSD(tgsd.P(i,1), tgsd.P(i,2), tgsd.P(i,3), tgsd.P(i,4), fullfile(run.name, 'TGSD', tgsd.N{i}));
    end
end

%% Configuration file
%elem    = {'plumeHeight', 'mass', 'duration', 'eddyConst', 'diffCoef', 'ftt', 'lithDens', 'pumDens', 'colStep', 'partStep', 'plume_model', 'plume_ratio', 'plumeAlpha', 'plumeBeta'};
% Append wind to
conf.windVel    = wind.vel; 
conf.windTrop   = wind.trop; 

var     = {'plumeHeight', 'mass', 'duration', 'eddyConst', 'diffCoef', 'ftt', 'lithDens', 'pumDens', 'colStep', 'partStep', 'plumeAlpha', 'plumeBeta', 'windVel','windTrop'};
varF    = {'Plume height', 'Mass', 'Duration', 'Eddy constant', 'Diffusion coef', 'FTT', 'Lithic density', 'Pumice density', 'Column steps', 'Particle steps', 'Alpha', 'Beta', 'wind', 'trop'};

% Index used to define what variables are associated for a given run.ass
if run.mass == 0
    idxE         = ~strcmp(var, 'duration');
elseif run.mass == 1
    idxE         = ~strcmp(var, 'mass');
end

% Here, generate a first batch of possibilities by considering permutations
% of wind velociy and duration
conf.elem    = var(idxE);
conf.elemF   = varF(idxE);
elemS   = elemC;

% [elements, elementsN, elemC, conf] ...
%     = prepareInput(elements, elementsN, elemS, elemC, conf, 'conf', '.conf');

[~,~,~, conf] ...
    = prepareInput(elements, elementsN, elemS, elemC, conf, 'conf', '.conf');

% Calculate MER
MER = zeros(size(conf.P,1),1);
for i = 1:size(conf.P,1)
    MER(i,1) = get_MER_DB12(conf.P(i, strcmp(var(idxE),'plumeHeight'))-vent.z, ...
        conf.P(i, strcmp(var(idxE),'windVel')),...
        conf.P(i, strcmp(var(idxE),'windTrop')));
end

% Remove the wind vectors associated with conf file from the global elements

conf.elem  = conf.elem(1:end-2);
conf.elemF = conf.elemF(1:end-2);

% Calculate mass and duration as a function of run.mass
if run.mass == 0
    duration    = conf.P(:, strcmp(var(idxE),'mass')) ./ MER ./ 3600;
    conf.P      = [conf.P, [MER, duration]];
    conf.elem   = [conf.elem, [{'MER'}, {'duration'}]];
    conf.elemF  = [conf.elemF, [{'MER'}, {'Duration'}]];
elseif run.mass == 1
    mass        = MER .* conf.P(:, strcmp(var(idxE),'windVel')).*3600;
    conf.P      = [conf.P, [MER, mass]];
    conf.elem   = [conf.elem, [{'MER'}, {'mass'}]];
    conf.elemF  = [conf.elemF, [{'MER'}, {'Mass'}]];
end
    
% generate a second bacth of possibilities for permutations on variables
% only relevant for the generation of configuration files
idx = ~(strcmp(conf.elem, 'MER') | strcmp(conf.elem, 'duration'));
conf.elem  = conf.elem(idx);
conf.elemF = conf.elemF(idx);

[elements, elementsN, elemC, conf] ...
    = prepareInput(elements, elementsN, elemS, elemC, conf, 'conf', '.conf');


% Stopped here. elements now contain all variable required to generate the
% permutations for conf files. Now still need to figure out how the indices
% of mass/height/wind can be used to later reshape and include MER and
% duration





















%%
elementsN = elementsN(~cellfun(@isempty, elements));
elements  = elements(~cellfun(@isempty, elements));
projectRun = makePossibilities(elements);
projectRun = mat2cell(projectRun, 192, ones(size(projectRun,2),1));



for i = 1:size(conf.P,1)
    makeConf(conf.P(i,1), conf.P(i,2), conf.P(i,3), tgsd.P(i,4), fullfile(run.name, 'TGSD', tgsd.N{i}));
end


% In case run mode is either single or sensitivity, remove reference to
% either mass or volume according to run.mass
% NEED TO ADAPT FOR PROBABILISTIC
if run.mode == 0 && run.mass == 1
    elem(2)     = [];
    elemF(2)    = [];
elseif run.mode == 0 && run.mass == 0  
    elem(3)     = [];
    elemF(3)    = []; 
end

for i = 1:length(elem)    
    elements{elemC}     = conf.(elem{i});
    elementsN{elemC}    = elemF{i};
    elemC               = elemC+1;
end

% Cleanup
elements  = elements(~cellfun(@isempty, elements));
elementsN = elementsN(~cellfun(@isempty, elements));




if run.mass == 0
    
end

% For each of grid, conf, gs and wind, return the all the combinations of
% variables and generate file names
function  [elements, elementsN, elemC, struc] = prepareInput(elements, elementsN, elemS, elemC, struc, fl, ext)
for i = 1:length(struc.elem)
    elements{elemC}     = struc.(struc.elem{i});
    elementsN{elemC}    = struc.elemF{i};
    elemC               = elemC+1;
end
struc.P  = makePossibilities(elements(elemS:elemC-1));   
struc.N  = cellfun(@num2str, num2cell(struc.P), 'UniformOutput', false);
struc.N  = makeName(struc.N, struc.elem, fl, ext);

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
combinations        = cellfun(@(x) x(:), combinations,'uniformoutput',false); %there may be a better way to do this
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


function makeConf

fid = fopen(fullfile(out_pth, 'CONF', seas_str{seas}, num2str(i), [num2str(j,'%04d'), '.conf']), 'wt');
fprintf(fid,...
    'PLUME_HEIGHT\t%d\nERUPTION_MASS\t%d\nVENT_EASTING\t%d\nVENT_NORTHING\t%d\nVENT_ELEVATION\t%d\nEDDY_CONST\t%d\nDIFFUSION_COEFFICIENT\t%d\nFALL_TIME_THRESHOLD\t%d\nLITHIC_DENSITY\t%d\nPUMICE_DENSITY\t%d\nCOL_STEPS\t%d\nPART_STEPS\t%d\nPLUME_MODEL\t%d\nALPHA\t%d\nBETA\t%d\n',...
    ht_tmp(j), mass_tmp(j),...
    data.vent_easting, data.vent_northing, data.vent_ht,...
    data.eddy_const, data.diff_coeff, data.ft_thresh,...
    data.lithic_dens, data.pumice_dens,...
    data.col_step, data.part_step, 2, data.alpha, data.beta);
fclose(fid);
% 
% 
% 
% 
% grid = [1000,5000];%,10000,20000];
% wind = 10:10:40;
% height = (10:5:30).*1000;
% MER = zeros(length(height), length(wind));
% T = 1*3600;
% 
% parfor iG = 1:length(grid)
%     for iH = 1:length(height)
%         for iW = 1:length(wind)
%                         
%             
%             gsFl    = 'GS_T2_D.gsd';
%             windFL  = ['wind_', num2str(wind(iW), '%2.0f'), '.wind'];
%             grdFL   = ['stromboli',num2str(grid(iG)),'.utm'];
% 
%             outFLt   = sprintf('out/TPg%3.0f_h%2.0f_w%2.0f.out', grid(iG), height(iH)./1e3, wind(iW) );
%             confFlt  = sprintf('conf/TPg%3.0f_h%2.0f_w%2.fd.conf', grid(iG), height(iH)./1e3, wind(iW) );
%             figFlt  = sprintf('fig/TPg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
%             
%             outFLg   = sprintf('out/GHg%3.0f_h%2.0f_w%2.0f.out', grid(iG), height(iH)./1e3, wind(iW) );
%             confFlg  = sprintf('conf/GHg%3.0f_h%2.0f_w%2.fd.conf', grid(iG), height(iH)./1e3, wind(iW) );
%             figFlg  = sprintf('fig/GHg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
%             
%             figFld  = sprintf('fig/DIFFg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
%             
%             conf = fileread('dummy.conf');
%             conf = strrep(conf, 'PLMHT', num2str(height(iH)));
%             conf = strrep(conf, 'MSS', num2str(get_mer(height(iH), wind(iW)) .* T ));
%             fid = fopen(confFlt, 'w');
%             fprintf(fid, '%s', conf);
%             fclose(fid);
%             
%             conf = fileread('dummyGH.conf');
%             conf = strrep(conf, 'PLMHT', num2str(height(iH)));
%             conf = strrep(conf, 'MSS', num2str(get_mer(height(iH), wind(iW)) .* T ));
%             fid = fopen(confFlg, 'w');
%             fprintf(fid, '%s', conf);
%             fclose(fid);
%             
%             alt  = 0:500:35000;
%             wnd = ones(length(alt),3);
%             wnd(:,1) = alt;
%             wnd(1,1) = 1;
%             wnd(:,3) = wnd(:,3).*90; 
%             wnd(1,2) = 1;
%             wnd(wnd(:,1)==10000,2) = wind(iW);
%             wnd(wnd(:,1)==20000,2) = 1;
%             wnd(:,2) = interp1([1,10000,20000,35000], [1, wind(iW), 1,1], wnd(:,1));
%             dlmwrite(windFL, wnd, 'delimiter', '\t', 'precision', 5);
%             
%             
%             modSt = sprintf('./tephra2-2012 %s %s %s %s > %s', confFlt, grdFL, windFL, gsFl, outFLt);
%             system(modSt);
%             modSg = sprintf('./tephra2-2012-GH %s %s %s > %s', confFlg, grdFL, windFL, outFLg);
%             system(modSg);
%             
%             
%             [X,Y,Mt] = plotT2(outFLt); title('TephraProb')
%             f1 = gcf; caxis([1e-3 4]);
%             [X,Y,Mg] = plotT2(outFLg); title('GitHub')
%             f2 = gcf; caxis([1e-3 4]);
% 
%             saveas(f1, figFlt);
%             saveas(f2, figFlg);
%             
%             f3 = figure;
%             pcolor(X,Y,Mg-Mt), shading flat
%             colorbar, title('Difference (kg m-2');
%             saveas(f3, figFld);
%             
%             close(f1)
%             close(f2)
%             close(f3)
%              % MER(iH,iW) = get_mer(height(iH), wind(iW));
%         end
%     end
% end
