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
force = 1; % Forces the code to re-generate individual figures
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
bashF      = prepareASCII(fullfile(fold, filesep, 'inversionConfig.conf'));
inversion.inFile.input = parseValue(bashF, 'inputFile=', 0);
inversion.inFile.wind  = parseValue(bashF, 'windFile=', 0);
inversion.inFile.grid  = parseValue(bashF, 'gridFile=', 0);
inversion.inFile.conf  = parseValue(bashF, 'configFile=', 0);
inversion.batch        = logical(str2double(parseValue(bashF, 'BATCH=', 0)));
inversion.wind          = str2double(parseValue(bashF, 'fixedWind=',0));
inversion.vent.easting  = str2double(parseValue(bashF, 'ventE=',0));
inversion.vent.northing = str2double(parseValue(bashF, 'ventN=',0));
inversion.vent.zone     = str2double(parseValue(bashF, 'ventZ=',0));

% Read input points
inversion.observations = dlmread([fold, filesep, inversion.inFile.input]);

for i = 1:length(all)
    if exist([fold, filesep, all(i).name, filesep, 'parameters.README'], 'file')  
        % Read output 
        outF    = prepareASCII([fold, filesep, all(i).name, filesep, 'parameters.README']);
        idx     = find(strcmp(outF, 'Parameter Ranges:'));
        inversion.fit(i).folder = all(i).name;
        
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
            plotT2(fullfile(fold, all(i).name, 'tephra2.out'), 'zone', inversion.vent.zone, 'vent', [inversion.vent.easting, inversion.vent.northing], 'plot', 'linear', 'points', inversion.observations(:,[1,2,4])); % -> add zone somewhere        
            saveFig(gcf, fullfile(fold, all(i).name, ['tephra2.', figFormat]), figFormat);
        end

        % Beta distribution of the plume
        if ~exist(fullfile(fold, all(i).name, ['plume.', figFormat]), 'file') || force == 1
            plotBetaPlume(inversion.fit(i).alpha, inversion.fit(i).beta, inversion.fit(i).plumeHeight)     
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
fitResult          = [[inversion.fit.fit]', [inversion.fit.plumeHeight]', [inversion.fit.mass]', [inversion.fit.alpha]', [inversion.fit.beta]', [inversion.fit.diffCoef]', [inversion.fit.ftt]', [inversion.fit.medPhi]', [inversion.fit.stdPhi]']; 
[fitResult,fitI]   = sortrows(fitResult);
fitResult          = [[1:numel(fitI)]', fitResult];
fitTable           = array2table(fitResult,'VariableNames', {'Nb', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi'});
inversion.fitTable = fitTable;
% Write it to an asci file
fid = fopen([fold, filesep, 'summary.txt'], 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Nb', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'MdPhi', 'SigPhi');
for i = 1:size(fitResult, 1)
    fprintf(fid, '%i\t%5.2f\t%5.0f\t%g\t%5.2f\t%5.2f\t%5.0f\t%5.0f\t%5.2f\t%5.2f\n', fitResult(i,1), fitResult(i,2), fitResult(i,3), fitResult(i,4), fitResult(i,5), fitResult(i,6), fitResult(i,7), fitResult(i,8), fitResult(i,9), fitResult(i,10));
end
fclose(fid);

% If batch, plot space
if size(fitResult,1) > 1
    plot_rmse(fitTable, fold, 0);
    
    refine = 0;
    if size(fitTable,1) > 10
        choice = questdlg('Refine?', ...
            'Refine', ...
            'Yes, thanks for asking','Do I look like this type of guy?','Yes, thanks for asking');

        switch choice
            case 'Yes, thanks for asking'
                refine = 1;
            case 'Do I look like this type of guy?'
                refine = 0;
        end
    end
    
    % Function to refine the inversion plots
    if refine == 1
        fprintf('\tRefine the inversion plot by discarding values above a given fit threshold\n')
        check = 0;
        while check == 0
            figure; 
            plot(log10(fitTable.Fit), 1:numel(fitTable.Fit), '-r.', 'LineWidth', 2,'MarkerEdgeColor','k');
            xlabel('Fit'); ylabel('Frequency'); 
            title('Click where to split the data. You must choose at least 10 points.');
            fitT        = ginput(1);
            [~,bisIdx]  = min(abs(log10(fitTable.Fit)-fitT(1)));
            close(gcf);
            if bisIdx<10
                errordlg('You must chose at least 10 points. Try again!')
            else
                check = 1;
            end
        end
        plot_rmse(fitTable(1:bisIdx,:), fold, 1);
    end
end




%% Write excel file
%data = sortrows(stor,2);
% Write folder
xlwrite([rootFold, filesep, flName], {fold}, sheetName, 'A1');          % Write folder
xlwrite([rootFold, filesep, flName], {datestr(now)}, sheetName, 'A2');  % Write date
xlwrite([rootFold, filesep, flName], {'Input file:', inversion.inFile.input}, sheetName, 'A3');  
xlwrite([rootFold, filesep, flName], {'Configuration file:', inversion.inFile.conf}, sheetName, 'A4');
xlwrite([rootFold, filesep, flName], {'Grid file:', inversion.inFile.grid}, sheetName, 'A5');
xlwrite([rootFold, filesep, flName], {'Fixed wind:', num2str(inversion.wind)}, sheetName, 'A6');
xlwrite([rootFold, filesep, flName], {'Wind file:', inversion.inFile.wind}, sheetName, 'A7');
%if exist('grid_file', 'var');    xlwrite([rootFold, filesep, flName], {'Grid file:', inversion.inFile.grid}, sheetName, 'A5'); end;
%if exist('wind_file', 'var');    xlwrite([rootFold, filesep, flName], {'Wind file:', inversion.inFile.wind}, sheetName, 'A7'); end;

% Write ranges
var_names = {'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'Median', 'Sigma'};

%writetable(ranges, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'B8');
xlwrite([rootFold, filesep, flName], var_names, sheetName, 'B8');
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
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height.png", "Height")']}, sheetName, 'B11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass.png", "Mass")']}, sheetName, 'C11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff.png", "Diff")']}, sheetName, 'F11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt.png", "FTT")']}, sheetName, 'G11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med.png", "Median")']}, sheetName, 'H11');

    if  inversion.batch && refine == 1
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height_fine.png", "Height")']}, sheetName, 'B12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass_fine.png", "Mass")']}, sheetName, 'C12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff_fine.png", "Diff")']}, sheetName, 'F12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt_fine.png", "FTT")']}, sheetName, 'G12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med_fine.png", "Median")']}, sheetName, 'H12');
    end


% Write fits
%writetable(tb, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'A14');
xlwrite([rootFold, filesep, flName], {'Fit' 'Height' 'Mass' 'Alpha' 'Beta' 'Diff' 'FTT' 'MdPhi' 'SigPhi'}, sheetName, 'A14');
xlwrite([rootFold, filesep, flName], fitResult(:,2:end), sheetName, 'A15');

% Write figure links
map_link = cell(size(fitResult,1), 3);
for i = 1:size(fitResult,1)
    rootdir = [fold, filesep, 'res_', num2str(log10(inversion.ranges(i).mass(1))), '_', num2str(inversion.ranges(i).plumeHeight(1)/1000), filesep];
    map_link{i,1} = ['=HYPERLINK("', rootdir, 'tephra2.pdf"; "Map")'];
    map_link{i,2} = ['=HYPERLINK("', rootdir, 'in_vs_out.png"; "In vs out")'];
    map_link{i,3} = ['=HYPERLINK("', rootdir, 'plume.png"; "Plume")'];
end
xlwrite([rootFold, filesep, flName], map_link, sheetName, 'J15');



function plot_rmse(data, fold, type)

data.Mass = log10(data.Mass);
data.Fit  = log10(data.Fit);

%% 1 Mass
h1 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1); plot_inversion(data.Mass, data.Height, data.Fit, 'log10 Mass (kg)', 'Plume height (km asl)');
subplot(2,2,2); plot_inversion(data.Mass, data.MdPhi, data.Fit, 'log10 Mass (kg)', 'Median \phi');
subplot(2,2,3); plot_inversion(data.Mass, data.FTT, data.Fit, 'log10 Mass (kg)', 'Ftt (s)');
subplot(2,2,4); plot_inversion(data.Mass, data.Diff, data.Fit, 'log10 Mass (kg)', 'Diffusion');

%% 3 Med
h3 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1); plot_inversion(data.MdPhi, data.Mass, data.Fit, 'Median \phi', 'log10 Mass (kg)');
subplot(2,2,2); plot_inversion(data.MdPhi, data.Height, data.Fit, 'Median \phi', 'Plume height (km asl)');
subplot(2,2,3); plot_inversion(data.MdPhi, data.FTT, data.Fit, 'Median \phi', 'Ftt (s)');
subplot(2,2,4); plot_inversion(data.MdPhi, data.Diff, data.Fit, 'Median \phi', 'Diffusion');

%% 4 FTT
h4 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1); plot_inversion(data.FTT, data.Mass, data.Fit, 'Ftt (s)', 'log10 Mass (kg)');
subplot(2,2,2); plot_inversion(data.FTT, data.Height, data.Fit, 'Ftt (s)', 'Plume height (km asl)');
subplot(2,2,3); plot_inversion(data.FTT, data.MdPhi, data.Fit, 'Ftt (s)', 'Median \phi');
subplot(2,2,4); plot_inversion(data.FTT, data.Diff, data.Fit, 'Ftt (s)', 'Diffusion');
%%
%% Save figures
if type == 0
    saveas(h1, [fold, filesep, 'mass.png']);
    %saveas(h2, [fold, filesep, 'height.png']);
    saveas(h3, [fold, filesep, 'med.png']);
    saveas(h4, [fold, filesep, 'ftt.png']);
    %saveas(h5, [fold, filesep, 'diff.png']);
else
    saveas(h1, [fold, filesep, 'mass_fine.png']);
    %saveas(h2, [fold, filesep, 'height_fine.png']);
    saveas(h3, [fold, filesep, 'med_fine.png']);
    saveas(h4, [fold, filesep, 'ftt_fine.png']);
    %saveas(h5, [fold, filesep, 'diff_fine.png']);
end
    
function plot_inversion(xdata, ydata, zdata, xlab, ylab)

if numel(xdata<10)
    maxData = numel(xdata);
else
    maxData = 10;
end

cont    = 15;
x_step  = (max(xdata) - min(xdata))/cont; 
y_step  = (max(ydata) - min(ydata))/cont;

x       = min(xdata):x_step:max(xdata);
y       = min(ydata):y_step:max(ydata);

[xi,yi] = meshgrid(x,y) ;
zi      = griddata(xdata,ydata,zdata,xi,yi) ;

surfc(xi(1,:),yi(:,1),zi);  
axis square, hold on, box on

z_norm = linspace(40, 10, length(xdata));
scatter3(xdata, ydata, zdata, z_norm, 'k', 'fill', 'o', 'MarkerEdgeColor', 'k')
text(xdata(1:maxData),ydata(1:maxData),zdata(1:maxData), cellstr(sprintfc('%d',[1:maxData]')), 'FontSize', 11, 'Color', 'r', 'FontWeight', 'bold');


xlabel(xlab);
ylabel(ylab);

c = colorbar;
set(get(c,'ylabel'),'string','Log10 RMSE','fontsize',10);

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

figure, hold on;
line([0,maxVal], [0,maxVal], 'Color', 'k', 'LineWidth',1);
scatter(sqrt(obs), sqrt(comp), 35, log10(sqrt(abs(obs-comp))), 'filled', 'MarkerEdgeColor', 'k');
axis square tight, box on

xlabel('Sqrt observed accumulation (kg/m^2)');
ylabel('Sqrt computed accumulation (kg/m^2)');

% Plot wind
function plotWind(height, dir, speed)
figure;
subplot(1,2,1)
plot(speed, height, '-r+', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
xlabel('Wind speed (m/s)');
ylabel('Height (m asl)')
a = subplot(1,2,2);
plot(dir, height, '-r+', 'LineWidth', 1, 'MarkerEdgeColor', 'k');
xlabel('Direction (degree from N)');
xlim([0,360])
a.XTick = [0,90,180,270];

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
    