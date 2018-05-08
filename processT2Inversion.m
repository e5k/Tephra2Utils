function processT2Inversion

% Display results from TEPHRA2 inversion done on a Baobab cluster
% Written by S. Biass, 2015
% Updates:
% Feb 2016: Added a dependency to xlwrite to solve the problem of writing
%           excel on mac and linux

% Add dependancies
addpath('Inversion/PostProcessing/xlwrite/');
javaaddpath('Inversion/PostProcessing/xlwrite/poi_library/poi-3.8-20120326.jar');
javaaddpath('Inversion/PostProcessing/xlwrite/poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('Inversion/PostProcessing/xlwrite/poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('Inversion/PostProcessing/xlwrite/poi_library/xmlbeans-2.3.0.jar');
javaaddpath('Inversion/PostProcessing/xlwrite/poi_library/dom4j-1.6.1.jar');
   
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

read_check  = 0;
all         = dir(fold);        % Batch folder
all         = all([all.isdir]); % Retrieve dir only
all         = all(3:end);       % Remove system dir
stor        = []; %zeros(length(all));
count = 0;
for i = 1:length(all)
    if all(i).isdir
        count = count+1;
        % Read conf file
        if read_check == 0
            fid = fopen([fold, filesep, all(i).name, filesep, 'tmp.conf']);
            tline = fgetl(fid);
            while ischar(tline)
                if ~isempty(regexp(tline, 'FIXED_WIND', 'once'))
                    wind = get_num(tline);
                end
                tline = fgetl(fid); 
            end
            fclose(fid);
        end
        
        % Read run.sh
        if read_check == 0
            fid = fopen([fold, filesep, 'runMass.sh']);
            tline = fgetl(fid);
            while ischar(tline)
                if ~isempty(regexp(tline, 'INPUT=', 'once'))
                    [~, idx] = regexp(tline, '=');
                    input_file = tline(idx+1:end);
                elseif ~isempty(regexp(tline, 'WIND=', 'once'))
                    [~, idx] = regexp(tline, '=');
                    wind_file = tline(idx+1:end);
                elseif ~isempty(regexp(tline, 'GRID=', 'once'))
                    [~, idx] = regexp(tline, '=');
                    grid_file = tline(idx+1:end);
                elseif ~isempty(regexp(tline, 'configFile=', 'once'))
                    [~, idx] = regexp(tline, '=');
                    conf_file = tline(idx+1:end);
                end
                tline = fgetl(fid); 
            end
            fclose(fid); 
        end
        read_check = 1;
        
        % Read output
        stor_line = zeros(1,38);
        fid = fopen([fold, filesep, all(i).name, filesep, 'parameters.README']);
        
        if fid ~= -1
            tline = fgetl(fid);
            while ischar(tline)
                if ~isempty(regexp(tline, 'FIT =', 'once'))
                    stor_line(2) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Modeled Values:', 'once'))
                    ctrl = 0;
                elseif ~isempty(regexp(tline, 'Parameter Ranges:', 'once'))
                    ctrl = 1;
                elseif ~isempty(regexp(tline, 'Max Column Height:', 'once')) && ctrl == 0;
                    stor_line(3) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Total Mass Ejected:', 'once')) && ctrl == 0;
                    stor_line(4) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Alpha Param:', 'once')) && ctrl == 0;
                    stor_line(5) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Beta Param:', 'once')) && ctrl == 0;
                    stor_line(6) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Diffusion Coefficient:', 'once')) && ctrl == 0;
                    stor_line(7) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Fall Time Threshold:', 'once')) && ctrl == 0;
                    stor_line(8) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Eddy Constant:', 'once')) && ctrl == 0;
                    stor_line(9) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Max Particle Size', 'once')) && ctrl == 0;
                    stor_line(10) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Min Particle Size:', 'once')) && ctrl == 0;
                    stor_line(11) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Median Size:', 'once')) && ctrl == 0;
                    stor_line(12) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Std. Dev.:', 'once')) && ctrl == 0;
                    stor_line(13) = get_num(tline);
                elseif ~isempty(regexp(tline, '(Easting)', 'once')) && ctrl == 0;
                    stor_line(14) = get_num(tline);
                elseif ~isempty(regexp(tline, '(Northing)', 'once')) && ctrl == 0;
                    stor_line(15) = get_num(tline);
                elseif ~isempty(regexp(tline, 'Maximum Column Height:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(16) = t1;
                    stor_line(17) = t2;
                elseif ~isempty(regexp(tline, 'Alpha Param:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(18) = t1;
                    stor_line(19) = t2;
                elseif ~isempty(regexp(tline, 'Beta Param:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(20) = t1;
                    stor_line(21) = t2;
                elseif ~isempty(regexp(tline, 'Diffusion Coefficient:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(22) = t1;
                    stor_line(23) = t2;
                elseif ~isempty(regexp(tline, 'Eddy Constant:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(24) = t1;
                    stor_line(25) = t2;
                elseif ~isempty(regexp(tline, 'Total Mass Ejected:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(26) = t1;
                    stor_line(27) = t2;
                elseif ~isempty(regexp(tline, 'Median Grain Size:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(28) = t1;
                    stor_line(29) = t2;
                elseif ~isempty(regexp(tline, 'Std. Dev. in particle diameter:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(30) = t1;
                    stor_line(31) = t2;
                elseif ~isempty(regexp(tline, 'Fall Time Threshold:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(32) = t1;
                    stor_line(33) = t2;
                elseif ~isempty(regexp(tline, 'Wind Speed:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(34) = t1;
                    stor_line(35) = t2;
                elseif ~isempty(regexp(tline, 'Wind Direction:', 'once')) && ctrl == 1;
                    [t1, t2] = get_num(tline);
                    stor_line(36) = t1;
                    stor_line(37) = t2;
                end           
                tline = fgetl(fid); 
            end

            % Plot in vs out
            points = dlmread([fold, filesep, all(i).name, filesep, 'tephra.out'], '', 1, 0);
            if ~isnan(points(1,4))
                stor_line(38) = 0;  
            else
                stor_line(38) = 1;
            end
            stor_line(1) = count;
            stor = [stor; stor_line];
            fclose(fid);
        end
       
    end
end
dlmwrite([fold, filesep, 'stor_out.dat'], stor, 'delimiter', '\t');

stor_xls = sortrows([stor(:,2), stor(:,3), stor(:,4), stor(:,5), stor(:,6), stor(:,7), stor(:,8), stor(:,12), stor(:,13)],1);


fid = fopen([fold, filesep, 'summary.txt'], 'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Fit', 'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'Md Phi', 'Sig Phi');
for i = 1:size(stor_xls, 1)
    fprintf(fid, '%5.2f\t%5.2f\t%g\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\t%5.2f\n', stor_xls(i,1), stor_xls(i,2), stor_xls(i,3), stor_xls(i,4), stor_xls(i,5), stor_xls(i,6), stor_xls(i,7), stor_xls(i,8), stor_xls(i,9));
end
fclose(fid);

if count > 1
%    plot_rmse(stor, fold, 0);

    data = sortrows(stor,2);

    choice = questdlg('Refine?', ...
        'Refine', ...
        'Yes, thanks for asking','Do I look like this type of guy?','Yes, thanks for asking');

    switch choice
        case 'Yes, thanks for asking'
            refine = 1;
        case 'Do I look like this type of guy?'
            refine = 0;
    end

    if refine == 1
        h = figure; 
        %hist(log10(data(:,2)), 50);
        %hold on;
        plot( log10(data(:,2)), 1:length(data(:,2)), '-r', 'LineWidth', 2);
        xlabel('Error');
        ylabel('Frequency');
        title('Select a minimum of 10 points!')
        h2 = msgbox(sprintf('1. Choose data\n2. Right click\n3. Create variable\n4. Leave ans as a variable name\n5. Close figure'));
        fprintf('\t1. Choose data\n\t2. Right click\n\t3. Create variable\n\t4. Leave ans as a variable name\n\t5. Close figure')
        waitfor(h2);
        brush on
        waitfor(h)
        data = data(1:size(ans,1),:);
        plot_rmse(data, fold, 1);
        
    else
        plot_rmse(data, fold, 0);
    end
end

clc
tb = table(stor_xls(:,1), stor_xls(:,2), stor_xls(:,3), stor_xls(:,4), stor_xls(:,5), stor_xls(:,6),stor_xls(:,7),stor_xls(:,8),stor_xls(:,9),'VariableNames', {'Fit' 'Height' 'Mass' 'Alpha' 'Beta' 'Diff' 'FTT' 'MdPhi' 'SigPhi'})



%% Write excel file
data = sortrows(stor,2);
% Write folder
xlwrite([rootFold, filesep, flName], {fold}, sheetName, 'A1');          % Write folder
xlwrite([rootFold, filesep, flName], {datestr(now)}, sheetName, 'A2');  % Write date
xlwrite([rootFold, filesep, flName], {'Input:', input_file}, sheetName, 'A3');  
xlwrite([rootFold, filesep, flName], {'Conf:', conf_file}, sheetName, 'A4'); 
if exist('grid_file', 'var');    xlwrite([rootFold, filesep, flName], {'Grid:', grid_file}, sheetName, 'A5'); end;
xlwrite([rootFold, filesep, flName], {'Fixed wind:', num2str(wind)}, sheetName, 'A6');
if exist('wind_file', 'var');    xlwrite([rootFold, filesep, flName], {'Wind:', wind_file}, sheetName, 'A7'); end;

% Write ranges
var_names = {'Height', 'Mass', 'Alpha', 'Beta', 'Diff', 'FTT', 'Median', 'Sigma'};

%writetable(ranges, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'B8');
xlwrite([rootFold, filesep, flName], var_names, sheetName, 'B8');
xlwrite([rootFold, filesep, flName], {'Min'}, sheetName, 'A9');
xlwrite([rootFold, filesep, flName], {'Max'}, sheetName, 'A10');
xlwrite([rootFold, filesep, flName], [[min(stor(:,16)); max(stor(:,17))], [min(stor(:,26)); max(stor(:,27))], [min(stor(:,18)); max(stor(:,19))], [min(stor(:,20)); max(stor(:,21))], [min(stor(:,22)); max(stor(:,23))],[min(stor(:,32)); max(stor(:,33))],[min(stor(:,28)); max(stor(:,29))],[min(stor(:,30)); max(stor(:,31))]], sheetName, 'B9');

if count > 1
    % Link to space images
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height.png"; "Height")']}, sheetName, 'B11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass.png"; "Mass")']}, sheetName, 'C11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff.png"; "Diff")']}, sheetName, 'F11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt.png"; "FTT")']}, sheetName, 'G11');
    xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med.png"; "Median")']}, sheetName, 'H11');

    if refine == 1
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'height_fine.png"; "Height")']}, sheetName, 'B12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'mass_fine.png"; "Mass")']}, sheetName, 'C12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'diff_fine.png"; "Diff")']}, sheetName, 'F12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'ftt_fine.png"; "FTT")']}, sheetName, 'G12');
        xlwrite([rootFold, filesep, flName], {['=HYPERLINK("', fold, filesep, 'med_fine.png"; "Median")']}, sheetName, 'H12');
    end
end

% Write fits
%writetable(tb, [rootFold, filesep, flName], 'sheet', sheetName, 'Range', 'A14');
xlwrite([rootFold, filesep, flName], {'Fit' 'Height' 'Mass' 'Alpha' 'Beta' 'Diff' 'FTT' 'MdPhi' 'SigPhi'}, sheetName, 'A14');
xlwrite([rootFold, filesep, flName], [stor_xls(:,1), stor_xls(:,2), stor_xls(:,3), stor_xls(:,4), stor_xls(:,5), stor_xls(:,6),stor_xls(:,7),stor_xls(:,8),stor_xls(:,9)], sheetName, 'A15');

% Write figure links
map_link = cell(size(stor_xls,1), 3);
for i = 1:size(stor_xls,1)
    map_link{i,1} = ['=HYPERLINK("', fold, filesep, 'res_', num2str(log10(data(i,26))), '_', num2str(data(i,16)/1000), filesep, 'tephra2.pdf"; "Map")'];
    map_link{i,2} = ['=HYPERLINK("', fold, filesep, 'res_', num2str(log10(data(i,26))), '_', num2str(data(i,16)/1000), filesep, 'in_vs_out.png"; "In vs out")'];
    map_link{i,3} = ['=HYPERLINK("', fold, filesep, 'res_', num2str(log10(data(i,26))), '_', num2str(data(i,16)/1000), filesep, 'plume.png"; "Plume")'];
end
xlwrite([rootFold, filesep, flName], map_link, sheetName, 'J15');


function varargout = get_num(in_str)
tmp = regexp(in_str, '\s', 'split');
k = 1;
for i = 1:length(tmp)
    if ~isnan(str2double(tmp{i}))
        varargout{k} = str2double(tmp{i});
        k = k+1;
    end
end

function plot_rmse(data, fold, type)

%data = data(find(log10(data(:,4))>8.5 & log10(data(:,4))<10.1),:);

data = sortrows(data,2);

mass    = log10(data(:,4));
height  = data(:,3)./1000;
rms     = log10(data(:,2));
% a       = data(:,5);
% b       = data(:,6);
diff    = data(:,7);
ftt     = data(:,8);
med     = data(:,12);

%% 1 Mass
h1 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot_inversion(mass, height, rms, 'log10 Mass (kgm^-^3)', 'Plume height (km asl)');
subplot(2,2,2);
plot_inversion(mass, med, rms, 'log10 Mass (kgm^-^3)', 'Median \phi');
subplot(2,2,3);
plot_inversion(mass, ftt, rms, 'log10 Mass (kgm^-^3)', 'Ftt (s)');
subplot(2,2,4);
plot_inversion(mass, diff, rms, 'log10 Mass (kgm^-^3)', 'Diffusion');

%% 2 Height
h2 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot_inversion(height, mass, rms, 'Plume height (km asl)', 'log10 Mass (kgm^-^3)');
subplot(2,2,2);
plot_inversion(height, med, rms, 'Plume height (km asl)', 'Median \phi');
subplot(2,2,3);
plot_inversion(height, ftt, rms, 'Plume height (km asl)', 'Ftt (s)');
subplot(2,2,4);
plot_inversion(height, diff, rms, 'Plume height (km asl)', 'Diffusion');

%% 3 Med
h3 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot_inversion(med, mass, rms, 'Median \phi', 'log10 Mass (kgm^-^3)');
subplot(2,2,2);
plot_inversion(med, height, rms, 'Median \phi', 'Plume height (km asl)');
subplot(2,2,3);
plot_inversion(med, ftt, rms, 'Median \phi', 'Ftt (s)');
subplot(2,2,4);
plot_inversion(med, diff, rms, 'Median \phi', 'Diffusion');

%% 4 FTT
h4 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot_inversion(ftt, mass, rms, 'Ftt (s)', 'log10 Mass (kgm^-^3)');
subplot(2,2,2);
plot_inversion(ftt, height, rms, 'Ftt (s)', 'Plume height (km asl)');
subplot(2,2,3);
plot_inversion(ftt, med, rms, 'Ftt (s)', 'Median \phi');
subplot(2,2,4);
plot_inversion(ftt, diff, rms, 'Ftt (s)', 'Diffusion');

%% 5 Diffusion coefficient
h5 = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot_inversion(diff, mass, rms, 'Diffusion', 'log10 Mass (kgm^-^3)');
subplot(2,2,2);
plot_inversion(diff, height, rms, 'Diffusion', 'Plume height (km asl)');
subplot(2,2,3);
plot_inversion(diff, med, rms, 'Diffusion', 'Median \phi');
subplot(2,2,4);
plot_inversion(diff, ftt, rms, 'Diffusion', 'Ftt (s)');

%% Save figures
if type == 0
    saveas(h1, [fold, filesep, 'mass.png']);
    saveas(h2, [fold, filesep, 'height.png']);
    saveas(h3, [fold, filesep, 'med.png']);
    saveas(h4, [fold, filesep, 'ftt.png']);
    saveas(h5, [fold, filesep, 'diff.png']);
else
    saveas(h1, [fold, filesep, 'mass_fine.png']);
    saveas(h2, [fold, filesep, 'height_fine.png']);
    saveas(h3, [fold, filesep, 'med_fine.png']);
    saveas(h4, [fold, filesep, 'ftt_fine.png']);
    saveas(h5, [fold, filesep, 'diff_fine.png']);
end
    
function plot_inversion(xdata, ydata, zdata, xlab, ylab)

cont    = 20;

x_step  = (max(xdata) - min(xdata))/cont; 
y_step  = (max(ydata) - min(ydata))/cont;

x       = min(xdata):x_step:max(xdata);
y       = min(ydata):y_step:max(ydata);

[xi,yi] = meshgrid(x,y) ;

zi      = griddata(xdata,ydata,zdata,xi,yi,'linear') ;

surfc(xi(1,:),yi(:,1),zi);
axis square, hold on, box on

z_norm = linspace(40, 10, length(xdata));
scatter3(xdata, ydata, zdata, z_norm, 'k', 'fill', 'o', 'MarkerEdgeColor', 'k')
for i = 1:10
    text(xdata(i), ydata(i),zdata(i), num2str(i), 'FontSize', 8, 'Color', 'r', 'FontWeight', 'bold');
end

xlabel(xlab);
ylabel(ylab);

c = colorbar;
set(get(c,'ylabel'),'string','Log10 RMSE','fontsize',10);



