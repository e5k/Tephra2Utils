% 

function [massfit, massT2] = fitDeposit_new(fl,Cpl,varargin)

noPlot  = 0;
ctr     = [];

if nnz(strcmpi(varargin, '-noplot'))
    noPlot = 1;
end

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'contours')
        ctr = varargin{i+1};
        i = i+1;
    end
end

% Retrieves T2 as a matrix
if noPlot == 1
    [X,Y,M]     = plotT2(fl, '-noplot');
else
    [X,Y,M]     = plotT2(fl);
end

% Defines isomass values
if isempty(ctr)
    ctr = logspace(-1,3,10);
end

% Contour deposit
C   = contourc(X(1,:),Y(:,1),10.^M, ctr);

% Check if contours are produced
if isempty(C)
    %varargout = makeOut(nan, [], nargout);
    massfit = nan;
    massT2  = nan;
    warning('All contours are empty, returning Nan');
    return
end

% Retrieve polygons from contour cells
[Cx,Cy,Cz]  = C2xyz(C);
% Check if length of vectors are at least 5 to define a polygon
% Check that at least 2 points can be defined
if nnz(cellfun(@length, Cx)>=5) < 2 || length(unique(Cz)) == 1
    %varargout = makeOut(nan, [], nargout);  
    massfit = nan;
    massT2  = nan;
    warning('All contours are empty, returning Nan');
    return
else
    idxCtr = cellfun(@length, Cx)>=5;
end    

% Retrieve area of isomass
Ca = zeros(nnz(idxCtr),1);
count = 1;
for i = 1:length(Cz)
    if idxCtr(i) == 1
        Ca(count) = polyarea(Cx{i}, Cy{i})/1e6;
        count     = count+1;
    end
end

% Calculates the mass returned by Tephra2
massT2  = nansum(nansum(10.^M)).* (X(1,2)-X(1,1))^2;
massfit = tephraFits(sqrt(Ca), Cz(idxCtr),  {'exponential', 'powerlaw'}, 'dataType', 'isomass', 'segments', [1,4], 'C', Cpl);



%varargout = makeOut(massT2, massfit, nargout);


function argout = makeOut(mT2, mE, narg)
argout = cell(narg, 1);
for i=1:narg
    if i == 1;      argout{i} = mT2;
    elseif i == 2;  argout{i} = mE;
    end
end
