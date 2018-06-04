% Inversion
function flCell = prepareASCII(flPath)
    fid     = fopen(flPath);
    flCell  = textscan(fid, '%s', 'delimiter', '\n', 'MultipleDelimsAsOne',1);
    flCell  = flCell{1};
    fclose(fid);