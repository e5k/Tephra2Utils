
grid = [1000];%,10000,20000];
wind = 10:10:40;
height = (10:5:30).*1000;
MER = zeros(length(height), length(wind));
T = 1*3600;

parfor iG = 1:length(grid)
    for iH = 1:length(height)
        for iW = 1:length(wind)
                        
            
            gsFl    = 'GS_T2_D.gsd';
            windFL  = ['wind_', num2str(wind(iW), '%2.0f'), '.wind'];
            grdFL   = ['stromboli',num2str(grid(iG)),'.utm'];

            outFLt   = sprintf('out/TPg%3.0f_h%2.0f_w%2.0f.out', grid(iG), height(iH)./1e3, wind(iW) );
            confFlt  = sprintf('conf/TPg%3.0f_h%2.0f_w%2.fd.conf', grid(iG), height(iH)./1e3, wind(iW) );
            figFlt  = sprintf('fig/TPg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
            
            outFLg   = sprintf('out/GHg%3.0f_h%2.0f_w%2.0f.out', grid(iG), height(iH)./1e3, wind(iW) );
            confFlg  = sprintf('conf/GHg%3.0f_h%2.0f_w%2.fd.conf', grid(iG), height(iH)./1e3, wind(iW) );
            figFlg  = sprintf('fig/GHg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
            
            figFld  = sprintf('fig/DIFFg%3.0f_h%2.0f_w%2.fd.jpg', grid(iG), height(iH)./1e3, wind(iW) );
            
            conf = fileread('dummy.conf');
            conf = strrep(conf, 'PLMHT', num2str(height(iH)));
            conf = strrep(conf, 'MSS', num2str(get_mer(height(iH), wind(iW)) .* T ));
            fid = fopen(confFlt, 'w');
            fprintf(fid, '%s', conf);
            fclose(fid);
            
            conf = fileread('dummyGH.conf');
            conf = strrep(conf, 'PLMHT', num2str(height(iH)));
            conf = strrep(conf, 'MSS', num2str(get_mer(height(iH), wind(iW)) .* T ));
            fid = fopen(confFlg, 'w');
            fprintf(fid, '%s', conf);
            fclose(fid);
            
            alt  = 0:500:35000;
            wnd = ones(length(alt),3);
            wnd(:,1) = alt;
            wnd(1,1) = 1;
            wnd(:,3) = wnd(:,3).*90; 
            wnd(1,2) = 1;
            wnd(wnd(:,1)==10000,2) = wind(iW);
            wnd(wnd(:,1)==20000,2) = 1;
            wnd(:,2) = interp1([1,10000,20000,35000], [1, wind(iW), 1,1], wnd(:,1));
            dlmwrite(windFL, wnd, 'delimiter', '\t', 'precision', 5);
            
            
            modSt = sprintf('./tephra2-2012 %s %s %s %s > %s', confFlt, grdFL, windFL, gsFl, outFLt);
            system(modSt);
            modSg = sprintf('./tephra2-2012-GH %s %s %s > %s', confFlg, grdFL, windFL, outFLg);
            system(modSg);
            
            
            [X,Y,Mt] = plotT2(outFLt); title('TephraProb')
            f1 = gcf; caxis([1e-3 4]);
            [X,Y,Mg] = plotT2(outFLg); title('GitHub')
            f2 = gcf; caxis([1e-3 4]);

            saveas(f1, figFlt);
            saveas(f2, figFlg);
            
            f3 = figure;
            pcolor(X,Y,Mg-Mt), shading flat
            colorbar, title('Difference (kg m-2');
            saveas(f3, figFld);
            
            close(f1)
            close(f2)
            close(f3)
             % MER(iH,iW) = get_mer(height(iH), wind(iW));
        end
    end
end
