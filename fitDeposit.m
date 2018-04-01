% 

function varargout = fitDeposit(fl,Cpl)

[X,Y,M]     = plotT2(fl, '-noplot');

%ctr = 10.^(linspace(min(min(M)), max(max(M)), 20));
%ctr = 10.^(-3:4);%linspace(min(min(10.^M)), max(max(10.^M)), 15);
% ctr = zeros(9,length(0:3));
% for i = 1:9
%     ctr(i,:) = i.*(10.^(0:3));
% end
% ctr = reshape(ctr,1,numel(ctr));

ctr = logspace(0,3,9);

C   = contourc(X(1,:),Y(:,1),10.^M, ctr);

% Check if contours are produced
if isempty(C)
    varargout = makeOut(nan, [], nargout);
    return
end

[Cx,Cy,Cz]  = C2xyz(C);
% Check if length of vectors are at least 5 to define a polygon
% Check that at least 2 points can be defined
if nnz(cellfun(@length, Cx)>=5) < 2 || length(unique(Cz)) == 1
    varargout = makeOut(nan, [], nargout);
    return
else
    idxCtr = cellfun(@length, Cx)>=5;
end    

Ca         = zeros(nnz(idxCtr),1);

count = 1;
for i = 1:length(Cz)
    if idxCtr(i) == 1
        Ca(count) = polyarea(Cx{i}, Cy{i})/1e6;
        count     = count+1;
    end
end

massT2  = nansum(nansum(10.^M)).* (X(1,2)-X(1,1))^2;

massfit = cell(length(Cpl),1);
for i = 1:numel(Cpl)
    massfit{i} = tephraFits(sqrt(Ca), Cz(idxCtr),  {'exponential', 'powerlaw'}, 'dataType', 'isomass', 'segments', [1,4], 'C', Cpl(i),'plotType','none');

end

varargout = makeOut(massT2, massfit, nargout);


function argout = makeOut(mT2, mE, narg)
argout = cell(narg, 1);
for i=1:narg
    if i == 1;      argout{i} = mT2;
    elseif i == 2;  argout{i} = mE;
    end
end
