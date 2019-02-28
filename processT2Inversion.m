function inversion = processT2Inversion(varargin)

% Display results from TEPHRA2 inversion done on a Baobab cluster
% Written by S. Biass, 2015
% Updates:
% Feb 2016: Added a dependency to xlwrite to solve the problem of writing
%           excel on mac and linux
% Jun 2018: Re-wrote most of it

% Configuration
%shFile = 'runInversion_PBS.sh'; % Name of the file used to start the inversion

% Optional parameters
force = 0; % Forces the code to re-generate individual figures
figFormat = 'pdf'; % Accepts 'pdf', 'eps', 'png', 'fig'

% Add dependencies
if isempty(regexp(matlabpath, 'xlwrite', 'ONCE'))
    addpath('dependencies/xlwrite/');
    javaaddpath('dependencies/xlwrite/poi_library/poi-3.8-20120326.jar');
    javaaddpath('dependencies/xlwrite/poi_library/poi-ooxml-3.8-20120326.jar');
    javaaddpath('dependencies/xlwrite/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
    javaaddpath('dependencies/xlwrite/poi_library/xmlbeans-2.3.0.jar');
    javaaddpath('dependencies/xlwrite/poi_library/dom4j-1.6.1.jar');
end
% Read input folder
fold        = uigetdir;
if fold == 0
    return
end

% Clean paths
foldpart    = strsplit(fold, filesep); 
sheetName   = foldpart(end);                            % Defines the sheet name
flName      = [foldpart{end-1}, '.xls'];                % Defines the name of the excel file
rootFold    = [filesep, fullfile(foldpart{2:end-2})];   % Root folder

% Define excel sheet
if exist([rootFold, filesep, flName], 'file') 
    [~, sheet] = xlsfinfo(fullfile(rootFold, flName));
    if length(unique(strcmp(sheet, sheetName))) > 1
        sheetName = [sheetName, '_', num2str(length(unique(strcmp(sheet, sheetName))))];
    end
end

%% Read inversion output
inversion   = struct;

% Retrieve path to all inversion folders
all         = dir(fold);        % Batch folder
all         = all([all.isdir]); % Retrieve dir only
all         = all(3:end);       % Remove system dir

% Read runMass.sh
bashF                   = prepareASCII(fullfile(fold, filesep, 'inversionConfig.conf'));
inversion.inFile.input  = parseValue(bashF, 'inputFile=', 0);
inversion.inFile.wind   = parseValue(bashF, 'windFile=', 0);
inversion.inFile.grid   = parseValue(bashF, 'gridFile=', 0);
inversion.inFile.conf   = parseValue(bashF, 'configFile=', 0);
inversion.batch         = str2double(parseValue(bashF, 'BATCH=', 0));
inversion.wind          = str2double(parseValue(bashF, 'fixedWind=',0));
inversion.maxWindDir    = str2double(parseValue(bashF, 'maxWindDir=',0));
inversion.minWindDir    = str2double(parseValue(bashF, 'minWindDir=',0));
inversion.vent.easting  = str2double(parseValue(bashF, 'ventE=',0));
inversion.vent.northing = str2double(parseValue(bashF, 'ventN=',0));
inversion.vent.zone     = str2double(parseValue(bashF, 'ventZ=',0));
inversion.vent.zone     = str2double(parseValue(bashF, 'ventZ=',0));

% Read input points
inversion.observations = dlmread([fold, filesep, inversion.inFile.input]);

for i = 1:length(all)
    fprintf('Reading %s\n', all(i).name)
    if i == 1
        % Read tmp.conf to retrieve the SEED value. Do that only once
        tmpF            = prepareASCII([fold, filesep, all(i).name, filesep, 'tmp.conf']);
        inversion.seed  = str2double(parseValue(tmpF, 'SEED ',0));
    end
    
    if exist([fold, filesep, all(i).name, filesep, 'parameters.README'], 'file')  
        % Read output 
        outF    = prepareASCII([fold, filesep, all(i).name, filesep, 'parameters.README']);
        idx     = find(strcmp(outF, 'Parameter Ranges:'));
        inversion.fit(i).folder = all(i).name;
        
        % Retrieve run number (important for batch = 2)
        if inversion.batch == 2
            tmp = strsplit(all(i).name, '_');
            inversion.fit(i).Nb = str2double(tmp{end});
        else
            inversion.fit(i).Nb = i;
        end
        
        % Retrieve fit values
        flVar = {'FIT =', 'Max Column Height:', 'Total Mass Ejected:', 'Alpha Param:', 'Beta Param:', 'Diffusion Coefficient:', 'Fall Time Threshold:', 'Eddy Constant:', 'Median Size:','Std. Dev.:'};
        stVar = {'fit', 'plumeHeight', 'mass', 'alpha', 'beta', 'diffCoef', 'ftt', 'eddy', 'medPhi', 'stdPhi'};
        for j = 1:length(flVar)
            inversion.fit(i).(stVar{j}) = parseValue(outF(1:idx),flVar{j},1); 
        end

        % Retrieve ranges
        flVar = {'Column Height:', 'Alpha Param:', 'Beta Param:', 'Diffusion Coefficient:', 'Eddy Constant:', 'Total Mass Ejected:', 'Median Grain Size:', 'Std. Dev. in particle diameter:', 'Fall Time Threshold:', 'Wind Speed:' 'Wind Direction:'};
        stVar = {'plumeHeight', 'alpha', 'beta', 'diffCoef', 'eddy', 'mass', 'medPhi', 'stdPhi', 'ftt', 'windSpeed', 'windDir'};
        for j = 1:length(flVar)
            inversion.ranges(i).(stVar{j}) = parseValue(outF(idx:end),flVar{j},1); 
        end

        % Read the model output
        tmpPoints = dlmread([fold, filesep, all(i).name, filesep, 'tephra.out'], '', 1, 0);
        inversion.fit(i).modelOut = tmpPoints(:,4);

        % Read output wind
        tmpWind = dlmread([fold, filesep, all(i).name, filesep, 'wind_levels.out'], '', 1, 0);
        inversion.fit(i).wind.height    = tmpWind(:,1);
        inversion.fit(i).wind.speed     = tmpWind(:,2);
        inversion.fit(i).wind.direction = tmpWind(:,3);

        %% Plot figures
        % Plot map
        if ~exist(fullfile(fold, all(i).name, ['tephra2.', figFormat]), 'file') || force == 1
            plotT2(fullfile(fold, all(i).name, 'tephra2.out'), 'vent', [inversion.vent.easting, inversion.vent.northing], 'plot', 'linear', 'points', inversion.observations(:,[1,2,4]), '-novis'); % -> add zone somewhere        
            saveFig(gcf, fullfile(fold, all(i).name, ['tephra2.', figFormat]), figFormat);
        end

        % Beta distribution of the plume
        if ~exist(fullfile(fold, all(i).name, ['plume.', figFormat]), 'file') || force == 1
            plotBetaPlume(inversion.fit(i).alpha, inversion.fit(i).beta, inversion.fit(i).plumeHeight, '-noplot')     
            saveFig(gcf, fullfile(fold, all(i).name, ['plume.', figFormat]), figFormat);
        end

        % Observed vs computed
        if ~exist(fullfile(fold, all(i).name, ['obsVScomp.', figFormat]), 'file') || force == 1
            obsVScomp(inversion.observations(:,4), inversion.fit(i).modelOut);     
            saveFig(gcf, fullfile(fold, all(i).name, ['obsVScomp.', figFormat]), figFormat);
        end

        % Wind
        if ~exist(fullfile(fold, all(i).name, ['wind.', figFormat]), 'file') || force == 1
            plotWind(inversion.fit(i).wind.height, inversion.fit(i).wind.direction, inversion.fit(i).wind.speed);     
            saveFig(gcf, fullfile(fold, all(i).name, ['wind.', figFormat]), figFormat);
        end
        
        % TGSD
        if ~exist(fullfile(fold, all(i).name, ['TGSD.', figFormat]), 'file') || force == 1
            plotTGSD(inversion.fit(i).medPhi, inversion.fit(i).stdPhi);     
            saveFig(gcf, fullfile(fold, all(i).name, ['TGSD.', figFormat]), figFormat);
        end

        % Clean folder
        fl2del = {'h_q_rmse1.dat', 'model.out', 'node_', 'plume.dat', 'plume2.dat', 'results', 'tephra2.out.xyz'};
        for j = 1:length(fl2del)
            if exist(fullfile(fold, all(i).name, fl2del{j}), 'file')
                delete(fullfile(fold, all(i).name, fl2del{j}));
            end
        end
    else
        rmdir([fold, filesep, all(i).name],'s');
    end
end

%% Prepare outputs
% Main inversion result table
fitResult          = [[inversion.fit.Nb]', [inversion.fit.fit]', [inversion.fit.plumeHeight]', [inversion.fit.mass]', [inversion.fit.alpha]', [inversion.fit.beta]', [inversion.fit.diffCoef]', [inversion.fit.ftt]', [inversion.fit.medPhi]', [inversion.fit.stdPhi]']; 
fitResult          = sortrows(fitResult,2);
%fitResult          = [[1:numel(fitI)]', fitResult];
fitTable           = array2table(fitResult,'VariableNames', {'Nb', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi'});
inversion.fitTable = fitTable;
% Write it to an asci file
fid = fopen([fold, filesep, 'summary.txt'], 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Nb', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi');
for i = 1:size(fitResult, 1)
    fprintf(fid, '%i\t%5.2f\t%5.0f\t%g\t%5.2f\t%5.2f\t%5.0f\t%5.0f\t%5.2f\t%5.2f\n', fitResult(i,1), fitResult(i,2), fitResult(i,3), fitResult(i,4), fitResult(i,5), fitResult(i,6), fitResult(i,7), fitResult(i,8), fitResult(i,9), fitResult(i,10));
end
fclose(fid);

save([fold, filesep, foldpart{end-1}, '_', foldpart{end}, '.mat'], 'inversion')
disp(inversion.fitTable);

% If batch, plot space
if inversion.batch == 1
    plotT2Inversion(inversion)
elseif inversion.batch == 2
    plotT2Inversion(inversion)
    
    % Plot matrix of scatter plots
    fitPM = fitResult;
    fitPM(:,4) = log10(fitPM(:,4));
    fitPM(:,3) = fitPM(:,3)./1e3;
    varName = {'Nb', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi'};
    
    figure, [~,ax]= plotmatrix(fitPM);
    for i = 1:size(ax,1)
        xlabel(ax(end,i), varName{i})
        ylabel(ax(i,1), varName{i})
    end
    saveas(gcf, [fold, filesep, foldpart{end-1}, '_', foldpart{end}, '.fig']);
end




%% Write excel file
%data = sortrows(stor,2);
% Write folder
xlwrite([rootFold, filesep, flName], {fold}, sheetName, 'A1');          % Write folder
xlwrite([rootFold, filesep, flName], {datestr(now)}, sheetName, 'A2');  % Write date
xlwrite([rootFold, filesep, flName], {'Seed:', inversion.seed}, sheetName, 'A3');
xlwrite([rootFold, filesep, flName], {'Fixed wind:', inversion.wind}, sheetName, 'A4');
if inversion.wind == 0
    xlwrite([rootFold, filesep, flName], {'Wind range:', inversion.minWindDir, inversion.maxWindDir}, sheetName, 'A5');
else
    xlwrite([rootFold, filesep, flName], {'Wind file:', inversion.inFile.wind}, sheetName, 'A5');
end

% Write ranges
var_names = {'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'Median', 'Sigma'};

%writetable(ranges, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'B8');
try
    xlwrite([rootFold, filesep, flName], var_names, sheetName, 'B8');
catch ME
    disp('There was an error writing the excel sheet. If you are updating an existing run, try and delete the existing sheet from the excel workbook first.')
    return
end

xlwrite([rootFold, filesep, flName], {'Min'}, sheetName, 'A9');
xlwrite([rootFold, filesep, flName], {'Max'}, sheetName, 'A10');
xlwrite([rootFold, filesep, flName], [...
    getBnd(cat(1,inversion.ranges.plumeHeight)), ...
    getBnd(cat(1,inversion.ranges.mass)), ...
    getBnd(cat(1,inversion.ranges.alpha)), ...
    getBnd(cat(1,inversion.ranges.beta)), ...
    getBnd(cat(1,inversion.ranges.diffCoef)), ...
    getBnd(cat(1,inversion.ranges.ftt)), ...
    getBnd(cat(1,inversion.ranges.medPhi)), ...
    getBnd(cat(1,inversion.ranges.stdPhi))], sheetName, 'B9');


    % Link to space images
%     xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height.png", "Height")']}, sheetName, 'B11');
%     xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass.png", "Mass")']}, sheetName, 'C11');
%     xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff.png", "Diff")']}, sheetName, 'F11');
%     xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt.png", "FTT")']}, sheetName, 'G11');
%     xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med.png", "Median")']}, sheetName, 'H11');
% 
%     if  inversion.batch && refine == 1
%         xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height_fine.png", "Height")']}, sheetName, 'B12');
%         xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass_fine.png", "Mass")']}, sheetName, 'C12');
%         xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff_fine.png", "Diff")']}, sheetName, 'F12');
%         xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt_fine.png", "FTT")']}, sheetName, 'G12');
%         xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med_fine.png", "Median")']}, sheetName, 'H12');
%     end


% Write fits
%writetable(tb, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'A14');
xlwrite([rootFold, filesep, flName], {'Run', 'Fit' 'Height' 'Mass' 'Alpha' 'Beta' 'Diff' 'FTT' 'MdPhi' 'SigPhi'}, sheetName, 'A14');
xlwrite([rootFold, filesep, flName], fitResult, sheetName, 'A15');

% Write figure links
map_link = cell(size(fitResult,1), 3);
for i = 1:size(fitResult,1)
%     rootdir = [fold, filesep, 'mass', num2str(log10(inversion.ranges(i).mass(1))), '_ht', num2str(inversion.ranges(i).plumeHeight(1)/1000), filesep];
%     map_link{i,1} = ['=HYPERLINK("file:', filesep, filesep, rootdir, 'tephra2.pdf", "Map")'];
%     map_link{i,2} = ['=HYPERLINK("file:', filesep, filesep, rootdir, 'obsVScomp.pdf", "In vs out")'];
%     map_link{i,3} = ['=HYPERLINK("file:', filesep, filesep, rootdir, 'plume.pdf", "Plume")'];
%     map_link{i,3} = ['=HYPERLINK("file:', filesep, filesep, rootdir, 'wind.pdf", "Wind")'];
end
xlwrite([rootFold, filesep, flName], map_link, sheetName, 'J15');





function R = getBnd(data)
    R = [min(data(:,1)); max(data(:,2))];
    
% Search values in input files
function out = parseValue(flCell, string, type)
% Type: 0 string
%       1 numeric
%       2 expect a range

[~, idxE]   = regexp(flCell,string,'once');
idxVal      = find(cellfun('isempty',idxE)==0);
if isempty(idxVal)
    out = [];
    return
end

if type == 0
    out = strtrim(flCell{idxVal}(idxE{idxVal}+1:end));
elseif type == 1
    out = cellfun(@str2double, regexp(flCell{idxVal}, '\s', 'split'));
    out = out(~isnan(out));
% elseif type == 2   
%     out = cellfun(@str2double, regexp(flCell{idxVal(2)}, '\s', 'split'));
%     out = out(~isnan(out));
end

% Plot observed vs computed
function obsVScomp(obs, comp)

maxVal = sqrt(max(max([obs,comp])));

figure('Visible', 'off'), hold on;
line([0,maxVal], [0,maxVal], 'Color', 'k', 'LineWidth',1);
scatter(sqrt(obs), sqrt(comp), 35, log10(sqrt(abs(obs-comp))), 'filled', 'MarkerEdgeColor', 'k');
axis square tight, box on

xlabel('Sqrt observed accumulation (kg/m^2)');
ylabel('Sqrt computed accumulation (kg/m^2)');

% Plot wind
function plotWind(height, dir, speed)
figure('Visible', 'off');
subplot(1,2,1)
plot(speed, height, '-r+', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
xlabel('Wind speed (m/s)');
ylabel('Height (m asl)')
a = subplot(1,2,2);
plot(dir, height, '-r+', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
xlabel('Direction (degree from N)');
xlim([0,360])
a.XTick = [0,90,180,270];

% Plot TGSD
function plotTGSD(mu, sigma)
range = -20:20;
dist  = normpdf(range, mu, sigma);
dist2 = normcdf(range, mu, sigma);
idx = dist>1e-3;
figure('Visible', 'off');
bar(range(idx), dist(idx).*100);
xlabel('\phi')
ylabel('Wt. %')
yyaxis right
plot(range(idx), dist2(idx).*100, '-k', 'Linewidth',1);


% Save figure
function saveFig(f, pth, figFormat)

if strcmp(figFormat, 'fig')
    saveas(f, pth);
else
    f.PaperPositionMode = 'auto';
    print(f, pth, ['-d', figFormat]);
end
close(f)

% Data mining
function flCell = prepareASCII(flPath)
    fid     = fopen(flPath);
    flCell  = textscan(fid, '%s', 'delimiter', '\n', 'MultipleDelimsAsOne',1);
    flCell  = flCell{1};
    fclose(fid);
     